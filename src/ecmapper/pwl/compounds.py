"""Compound identity and BRENDA reaction equations for the PWL generator.

The BRENDA map export (``BRENDA_map_reactions.json``) gives every drawn
metabolite a display name, a BRENDA name and a BRENDA ligand ID. None of them
identifies a compound alone: the same compound has several display names
("Ac-CoA", "acetyl-CoA") and sometimes several ligand IDs (O2: 9 and 1180),
and some ligand IDs are wrong (NADH drawn with the ligand ID of H+).
:class:`CompoundRegistry` merges ligand IDs and BRENDA names into compounds,
see its docstring.
"""

from __future__ import annotations

import difflib
import re
from collections import Counter, defaultdict
from dataclasses import dataclass, field

from ..brenda_map import _normalize, expand_names
from .syntax import KEYWORDS

_LOCATION = re.compile(r"\s*(?:\[(side \d+|in|out|inside|outside)\]|/(in|out))$", re.IGNORECASE)
_PROTONS = {"h+", "h", "proton", "hydron"}


def split_location(name: str) -> tuple[str, str | None]:
    """Name without compartment and the compartment token: "H+[side 2]" gives ("H+", "side2")."""
    m = _LOCATION.search(name or "")
    if not m:
        return name, None
    loc = (m.group(1) or m.group(2)).lower().replace(" ", "")
    loc = {"inside": "in", "outside": "out"}.get(loc, loc)
    return name[: m.start()], loc


def compartment_label(loc: str) -> str:
    m = re.fullmatch(r"side(\d+)", loc)
    return f"side {m.group(1)}" if m else loc


_RING_FORM = re.compile(r"(?<=[a-z])o(?:pyranose|furanose)", re.IGNORECASE)


def _base(name: str) -> str:
    # ring forms of sugars are named like the sugar: beta-D-fructofuranose 6-phosphate ~ D-fructose 6-phosphate
    return _RING_FORM.sub("ose", split_location(name)[0])


def strict_key(name: str) -> str:
    """Name key that keeps stereo descriptors (D-/L-, (S)-, alpha-/beta-), see :func:`ecmapper.brenda_map.name_key`."""
    return _normalize(_base(name), stereo=True)


def loose_key(name: str) -> str:
    """Name key without stereo descriptors (and without ring forms of sugars)."""
    return _normalize(_base(name))


def slug(name: str) -> str:
    """A PWL token for a name: a trailing charge '+' or '-' is dropped ("NAD+" -> NAD), other characters become '_'."""
    n = (name or "").strip()
    n = re.sub(r"(?<=[A-Za-z0-9)\]])[+-]$", "", n)
    n = n.replace("+", "plus")
    s = re.sub(r"[^A-Za-z0-9]+", "_", n).strip("_") or "X"
    if s in KEYWORDS or s == "_":
        s += "_c"
    return s


@dataclass
class CompoundInfo:
    """A compound of the registry."""

    cid: int
    token: str = ""
    ligand_ids: set[int] = field(default_factory=set)
    display: Counter = field(default_factory=Counter)
    brenda_names: Counter = field(default_factory=Counter)
    reassigned: Counter = field(default_factory=Counter)  # (ligand id, display name) not merged into the compound

    @property
    def label(self) -> str:
        return self.display.most_common(1)[0][0] if self.display else self.token

    @property
    def count(self) -> int:
        return sum(self.display.values())

    def names(self) -> set[str]:
        return set(self.display) | set(self.brenda_names)


def _entries(data: dict):
    for r in data["reactions"]:
        for m in r["main_metabolites"]:
            yield r, "main", m
        for role in ("consumed", "produced", "undirected"):
            for m in r["cofactors"].get(role, []):
                yield r, role, m


class CompoundRegistry:
    """The compounds of a BRENDA map export, with one PWL token per compound for the whole map.

    Ligand IDs and BRENDA names (without compartment, stereo descriptors kept)
    are the nodes of a graph. A ligand ID and a name are joined if more than
    half of the occurrences of the ligand ID carry this name or more than half
    of the occurrences of the name carry this ligand ID. The connected
    components are the compounds. A drawn metabolite belongs to the compound of
    its BRENDA name, so a wrong ligand ID of a few metabolites (NADH with the
    ligand ID of H+) does not merge two compounds; it is listed in
    ``reassigned``. Two ligand IDs of the same name (O2: 9 and 1180) end up in
    one compound.

    The token is made from the most frequent display name (:func:`slug`);
    compounds are numbered by frequency, so the more frequent one keeps the
    plain token if two would get the same one, the others get ``_L<ligand ID>``.
    """

    def __init__(self, data: dict, synonyms: dict[str, tuple[str, ...]] | None = None):
        parent: dict = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        def union(a, b):
            ra, rb = find(a), find(b)
            if ra != rb:
                parent[ra] = rb

        pair, by_lig, by_name = Counter(), Counter(), Counter()
        occurrences = []
        for _r, _role, m in _entries(data):
            display = split_location(m.get("display_name") or "")[0]
            bname = split_location(m.get("brenda_name") or m.get("display_name") or "")[0]
            key = ("N", strict_key(bname) or strict_key(display) or display)
            lig = m.get("brenda_ligand_id")
            find(key)
            occurrences.append((lig, key, display, bname))
            by_name[key] += 1
            if lig is not None:
                pair[(lig, key)] += 1
                by_lig[lig] += 1
        for (lig, key), n in pair.items():
            if n > 0.5 * by_lig[lig] or n > 0.5 * by_name[key]:
                union(("L", lig), key)
        roots: dict = {}
        self._info: list[CompoundInfo] = []
        self._by_key: dict = {}
        for lig, key, display, bname in occurrences:
            root = find(key)
            if root not in roots:
                roots[root] = len(self._info)
                self._info.append(CompoundInfo(len(self._info)))
            info = self._info[roots[root]]
            info.display[display] += 1
            info.brenda_names[bname] += 1
            if lig is not None:
                if find(("L", lig)) == root:
                    info.ligand_ids.add(lig)
                else:
                    info.reassigned[(lig, display)] += 1
            self._by_key[key] = info.cid
        self._assign_tokens()
        self.set_synonyms(synonyms)

    def set_synonyms(self, synonyms: dict[str, tuple[str, ...]] | None) -> None:
        """Use a synonym list {identifier: names}, e.g. :class:`~ecmapper.sources.KEGGCompoundNames` data.

        Names of compounds and of equation terms that are a name of the same
        identifier match; a compound whose names lead to exactly one identifier
        gets it as cross-reference (:meth:`synonym_ids`).
        """
        self._syn: dict[str, set[str]] = {}
        for identifier, names in (synonyms or {}).items():
            for name in names:
                self._syn.setdefault(strict_key(name), set()).add(identifier)
        self._keys = [self._compound_keys(info) for info in self._info]

    def synonym_ids(self, name: str) -> set[str]:
        return set(self._syn.get(strict_key(name), ()))

    def _assign_tokens(self) -> None:
        taken: set[str] = set()
        for info in sorted(self._info, key=lambda i: (-i.count, i.label)):
            base = slug(info.label)
            token = base
            if token in taken:
                token = f"{base}_L{min(info.ligand_ids)}" if info.ligand_ids else f"{base}_{info.cid}"
            k = 2
            while token in taken:
                token = f"{base}_{k}"
                k += 1
            taken.add(token)
            info.token = token
        self._by_token = {info.token: info for info in self._info}

    def _compound_keys(self, info: CompoundInfo) -> tuple[set[str], set[str], set[str]]:
        names = expand_names(info.names())
        ids = set().union(set(), *(self._syn.get(strict_key(n), set()) for n in names))
        return {strict_key(n) for n in names} - {""}, {loose_key(n) for n in names} - {""}, ids

    def xref_ids(self, info: CompoundInfo) -> set[str]:
        """Identifiers of the synonym list for a compound: all its names lead to the same single one."""
        ids = self._keys[info.cid][2]
        return ids if len(ids) == 1 else set()

    # --- lookup
    def compound(self, entry: dict) -> CompoundInfo:
        """The compound of a metabolite entry of the export (main metabolite or cofactor)."""
        display = split_location(entry.get("display_name") or "")[0]
        bname = split_location(entry.get("brenda_name") or entry.get("display_name") or "")[0]
        key = ("N", strict_key(bname) or strict_key(display) or display)
        return self._info[self._by_key[key]]

    def token(self, entry: dict) -> str:
        return self.compound(entry).token

    def by_token(self, token: str) -> CompoundInfo:
        return self._by_token[token]

    def keys(self, info: CompoundInfo) -> tuple[set[str], set[str], set[str]]:
        """(strict, loose) name keys of a compound and its identifiers in the synonym list."""
        return self._keys[info.cid]

    def __iter__(self):
        return iter(sorted(self._info, key=lambda i: i.token))

    def __len__(self) -> int:
        return len(self._info)

    def to_json(self) -> dict:
        return {
            info.token: {
                "label": info.label,
                "brenda_ligand_ids": sorted(info.ligand_ids),
                "display_names": dict(info.display.most_common()),
                "brenda_names": dict(info.brenda_names.most_common()),
                "reassigned": [{"ligand_id": lig, "display_name": d, "count": n}
                               for (lig, d), n in sorted(info.reassigned.items(), key=str)],
            }
            for info in self
        }


# ---------------------------------------------------------------------------------- equations

_COEF = re.compile(r"^(\d+(?:\.\d+)?|n)\s+(?=\S)")


@dataclass
class Term:
    name: str
    coef: str | None = None
    location: str | None = None


def parse_equation(equation: str) -> tuple[list[Term], list[Term]] | None:
    """Sides of a BRENDA equation ``2 A + B[side 1] = C``; None if it is not of this form."""
    if equation.count(" = ") != 1:
        return None
    sides = []
    for side in equation.split(" = "):
        terms = []
        for raw in side.split(" + "):
            raw = raw.strip()
            if not raw:
                return None
            coef = None
            m = _COEF.match(raw)
            if m:
                coef, raw = m.group(1), raw[m.end():]
            name, loc = split_location(raw)
            terms.append(Term(name.strip(), coef, loc))
        sides.append(terms)
    return sides[0], sides[1]


def _term_keys(term: Term) -> tuple[set[str], set[str]]:
    names = expand_names([term.name])
    return {strict_key(n) for n in names} - {""}, {loose_key(n) for n in names} - {""}


def is_proton(name: str) -> bool:
    return strict_key(name) in _PROTONS


@dataclass
class EquationMatch:
    """The sides of the participants of a drawn reaction in a BRENDA equation.

    `sides` maps the index of a participant to "L" or "R" (participants on
    both sides or not found are left out); `coef` and `location` hold the
    coefficient and compartment of the matching term.
    """

    equation: str
    reaction_id: str
    sides: dict[int, str]
    coef: dict[int, str]
    location: dict[int, str]
    score: int
    both: set[int] = field(default_factory=set)
    methods: dict[int, str] = field(default_factory=dict)  # participant -> "kegg", "single" or "similar"


def match_equation(registry: CompoundRegistry, participants: list[CompoundInfo], equation: str, reaction_id: str,
                   main: list[int] | None = None):
    """Match the participants (compounds of a drawn reaction) with the terms of an equation.

    `main` are the indices of the main metabolites among the participants; if
    exactly one of them and exactly one term of the equation are not matched,
    they are taken to be the same compound (BRENDA uses other names than the
    map for some compounds, e.g. "glycerone phosphate" for DHAP).
    """
    parsed = parse_equation(equation)
    if parsed is None:
        return None
    keys = [registry.keys(c) for c in participants]
    methods: dict[int, str] = {}
    sides: dict[int, set[str]] = defaultdict(set)
    coef: dict[int, str] = {}
    loc: dict[int, str] = {}
    term_index: dict[int, int] = {}
    for side, terms in zip("LR", parsed):
        for term in terms:
            term_index[id(term)] = 0
            strict, loose = _term_keys(term)
            hits = [i for i, (s, _l, _x) in enumerate(keys) if s & strict]
            if not hits:
                loose_hits = [i for i, (_s, lo, _x) in enumerate(keys) if lo & loose]
                # without stereo descriptors only if unambiguous among the participants
                if len({participants[i].cid for i in loose_hits}) == 1:
                    hits = loose_hits
            if not hits:
                ids = registry.synonym_ids(term.name)
                syn_hits = [i for i, (_s, _l, x) in enumerate(keys) if x & ids]
                if len({participants[i].cid for i in syn_hits}) == 1:
                    hits = syn_hits
                    methods.update({i: "kegg" for i in hits})
            term_index[id(term)] = len(hits)
            for i in hits:
                sides[i].add(side)
                if term.coef and term.coef != "1":
                    coef[i] = term.coef
                if term.location:
                    loc[i] = term.location
    if main is not None:
        # exactly one main metabolite and exactly one term (no proton) are left: they are the same compound
        free_main = [i for i in main if i not in sides]
        free_terms = [(side, term) for side, terms in zip("LR", parsed) for term in terms
                      if not is_proton(term.name) and term_index[id(term)] == 0]
        pairs, method = [], None
        if len(free_main) == 1 and len(free_terms) == 1:
            pairs, method = [(free_main[0], free_terms[0])], "single"
        elif free_main and len(free_terms) >= len(free_main):
            pairs, method = _similar_pairs(registry, participants, free_main, free_terms), "similar"
        for i, (side, term) in pairs:
            methods[i] = method
            sides[i].add(side)
            if term.coef and term.coef != "1":
                coef[i] = term.coef
            if term.location:
                loc[i] = term.location
    one = {i: next(iter(s)) for i, s in sides.items() if len(s) == 1}
    both = {i for i, s in sides.items() if len(s) == 2}
    return EquationMatch(equation, reaction_id, one, {i: c for i, c in coef.items() if i in one},
                         {i: x for i, x in loc.items()}, len(one), both, methods)


def _similarity(a: str, b: str) -> float:
    return difflib.SequenceMatcher(None, a, b).ratio()


def _similar_pairs(registry, participants, free_main, free_terms, minimum: float = 0.6, margin: float = 0.15):
    """Main metabolites and terms with similar names (by loose key): each the clearly best partner of the other."""
    score = {}
    for i in free_main:
        names = {loose_key(n) for n in participants[i].names()} - {""}
        for k, (_side, term) in enumerate(free_terms):
            score[(i, k)] = max((_similarity(n, loose_key(term.name)) for n in names), default=0.0)
    pairs = []
    for i in free_main:
        ranked = sorted(((score[(i, k)], k) for k in range(len(free_terms))), reverse=True)
        best, k = ranked[0]
        second = ranked[1][0] if len(ranked) > 1 else 0.0
        others = [score[(j, k)] for j in free_main if j != i]
        if best >= minimum and best - second >= margin and all(best - o >= margin for o in others):
            pairs.append((i, free_terms[k]))
    return pairs


def equation_signature(registry: CompoundRegistry, equation: str) -> tuple[frozenset, frozenset] | None:
    """The sides of an equation as sets of loose name keys (without protons), to compare directions."""
    parsed = parse_equation(equation)
    if parsed is None:
        return None
    return tuple(frozenset(loose_key(t.name) for t in side if not is_proton(t.name)) for side in parsed)
