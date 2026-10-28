"""Find the reactions of a model on a BRENDA pathway map."""

from __future__ import annotations

import html
import json
import os
import re
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

import pandas as pd

from . import annotation
from .disambiguation import covers
from .ec import is_complete, normalize_ec
from .mapper import ECMapper
from .metabolites import CURRENCY, PROTON, MetaboliteComparison
from .sources import (
    BiGGrMetabolites,
    BKMSSource,
    CachedResource,
    KEGGCompoundNames,
    KEGGReactionCompounds,
    KEGGSource,
)

#: values of :attr:`ReactionMapMatch.status`
MAP_STATUSES = ("on map", "related", "not on map", "boundary")

#: modes of the search in the vicinity of found nodes, see :meth:`MapMatcher.search_vicinity`
VICINITY_MODES = ("neighbours", "propagate")

_STEREO = re.compile(r"(?<![a-z0-9])(?:all-trans|\((?:[0-9]*[rsez][,0-9rsez]*)\)|dl|[dl]|alpha|beta|cis|trans)-")
_PROTON_NAMES = {"h+", "h", "proton", "protons", "hydron", "hydrogenion"}
# compartments in names of map metabolites, e.g. "H+[side 1]"
_SIDE = re.compile(r"\s*\[(?:side \d+|in|out|inside|outside)\]$")

#: generic classes of redox carriers {class: pattern of name keys}; maps often show "a quinone",
#: models specific ones like ubiquinone-8
GROUPS = {
    "quinone": re.compile(r"(?:a|an|\d*)(?:ubi|mena|demethylmena|plasto)?quinone\d*|coenzymeq\d*|oxidizedquinone"),
    "quinol": re.compile(
        r"(?:a|an|\d*)(?:ubi|mena|demethylmena|plasto)?quinol\d*|hydroquinone|reduced(?:ubi|mena)?quinone\d*"
    ),
}

#: keys of metabolites that count as currency metabolites in the score
CURRENCY_KEYS = CURRENCY | {f"group:{g}" for g in GROUPS}


def name_key(name: str) -> str:
    """Key for comparing metabolite names.

    Case, punctuation and stereo descriptors (``D-``, ``(S)-``, ``alpha-``,
    ...) are ignored and ``bisphosph...``/``diphosph...`` as well as ``... acid``/``...ate``
    are treated alike, so that e.g. "fructose-1,6-bisphosphate" and
    "D-Fructose 1,6-diphosphate" get the same key. D- and L-forms of a
    compound therefore get the same key, too. Sides of membranes ("H+[side 1]")
    are removed.
    """
    n = _SIDE.sub("", name.lower())
    n = _STEREO.sub("", n)
    n = n.replace("bisphosph", "diphosph").replace("ic acid", "ate")
    return re.sub(r"[^a-z0-9+]", "", n)


_EITHER = re.compile(r"\(([a-z])\)", re.IGNORECASE)


def expand_names(names: Iterable[str]) -> set[str]:
    """Names with optional parts in both forms: "NAD(P)H" gives NADH and NADPH as well."""
    result = set()
    for name in names:
        if not name:
            continue
        result.add(name)
        if _EITHER.search(name):
            result |= {_EITHER.sub("", name), _EITHER.sub(r"\1", name)}
    return result


def _keys(compounds: Iterable[str], names: Iterable[str]) -> frozenset[str]:
    """Keys of a metabolite: KEGG compound IDs, ``name:`` keys and ``group:`` keys of generic classes."""
    names = {n for n in names if n and n not in _PROTON_NAMES}
    groups = {f"group:{g}" for g, pattern in GROUPS.items() if any(pattern.fullmatch(n) for n in names)}
    return frozenset(compounds) | frozenset(f"name:{n}" for n in names) | frozenset(groups)


@dataclass
class MapMetabolite:
    """A metabolite drawn next to an enzyme node."""

    names: tuple[str, ...]
    ligand_id: int | None = None
    role: str = "main"  # main, consumed or produced
    node_id: str | None = None  # the drawn metabolite node, may be shared with other enzyme nodes


@dataclass
class MapNode:
    """An enzyme node of the map."""

    node_id: str
    subsystem: str
    ec_label: str
    ec_numbers: list[str]
    enzyme_name: str = ""
    brenda_reaction_ids: list[str] = field(default_factory=list)
    metabolites: list[MapMetabolite] = field(default_factory=list)
    x: float | None = None
    y: float | None = None


class BrendaMap:
    """The enzyme nodes of a BRENDA pathway map.

    `links` are pairs of drawn metabolite nodes connected by a line, e.g. the
    same metabolite in two subsystems.
    """

    def __init__(self, nodes: Iterable[MapNode], links: Iterable[tuple[str, str]] = ()):
        self.nodes = {n.node_id: n for n in nodes}
        self.links = list(links)
        self._places: dict[str, list[tuple[str, MapMetabolite]]] | None = None
        self._neighbours: dict[str, list[tuple[str, str]]] = {}

    def __len__(self) -> int:
        return len(self.nodes)

    def _build_places(self) -> None:
        drawn: dict[str, list[tuple[str, MapMetabolite]]] = {}
        for node in self.nodes.values():
            for m in node.metabolites:
                if m.role == "main" and m.node_id:
                    drawn.setdefault(m.node_id, []).append((node.node_id, m))
        parent: dict[str, str] = {}

        def root(x: str) -> str:
            while parent.get(x, x) != x:
                parent[x] = parent.get(parent[x], parent[x])
                x = parent[x]
            return x

        for a, b in self.links:
            if a in drawn and b in drawn and root(a) != root(b):
                parent[root(a)] = root(b)
        places: dict[str, list[tuple[str, MapMetabolite]]] = {}
        for node_id in sorted(drawn):
            places.setdefault(root(node_id), []).extend(drawn[node_id])
        neighbours: dict[str, list[tuple[str, str]]] = {}
        for place, items in places.items():
            enzymes = list(dict.fromkeys(e for e, _ in items))
            for a in enzymes:
                neighbours.setdefault(a, []).extend((b, place) for b in enzymes if b != a)
        self._places, self._neighbours = places, neighbours

    def neighbours(self, node_id: str) -> list[tuple[str, str]]:
        """The enzyme nodes next to a node, as (enzyme node, place) for each main metabolite they share.

        A place is a drawn main metabolite node together with the metabolite
        nodes linked to it; it is identified by one of their IDs.
        """
        if self._places is None:
            self._build_places()
        return self._neighbours.get(node_id, [])

    def place_metabolites(self, place: str) -> list[MapMetabolite]:
        """The metabolites of the enzyme nodes drawn at a place."""
        if self._places is None:
            self._build_places()
        return [m for _, m in self._places.get(place, [])]

    @classmethod
    def from_json(cls, path: str | os.PathLike) -> BrendaMap:
        """Read the reactions export of a map (``reactions`` with EC numbers, metabolites and layout)."""
        with open(path, encoding="utf-8") as f:
            data = json.load(f)
        nodes = []
        for r in data["reactions"]:
            metabolites = [
                MapMetabolite(
                    tuple(dict.fromkeys(n for n in (m.get("brenda_name"), m.get("display_name")) if n)),
                    m.get("brenda_ligand_id"),
                    role,
                    m.get("node_id"),
                )
                for role, items in (
                    ("main", r.get("main_metabolites", [])),
                    ("consumed", r.get("cofactors", {}).get("consumed", [])),
                    ("produced", r.get("cofactors", {}).get("produced", [])),
                    ("undirected", r.get("cofactors", {}).get("undirected", [])),
                )
                for m in items
            ]
            layout = r.get("layout", {})
            nodes.append(
                MapNode(
                    r["node_id"],
                    r.get("subsystem", ""),
                    r.get("ec_label", ""),
                    list(r.get("ec_numbers", [])),
                    r.get("enzyme_name", ""),
                    list(r.get("brenda_reaction_ids", [])),
                    metabolites,
                    layout.get("x"),
                    layout.get("y"),
                )
            )
        links = [
            (link["start_node_id"], link["end_node_id"])
            for link in data.get("metabolite_links", []) + data.get("unassigned_links", [])
            if link.get("start_node_id") and link.get("end_node_id")
        ]
        return cls(nodes, links)


@dataclass
class NodeMatch:
    """A map node that is a candidate for a model reaction.

    `via` lists how the node was found: ``reaction-id`` (a BRENDA reaction of
    the node is aligned in BKMS-react with a KEGG or MetaCyc reaction of the
    model reaction), ``ec`` (common EC number), ``metabolites`` (at least two
    common main metabolites) and ``vicinity`` (accepted because it is next to a
    node of a model reaction, see :meth:`MapMatcher.search_vicinity`; that
    match is `vicinity_of`, as (reaction ID, node ID)). `score` compares the
    metabolites of the model reaction with those of the node (main metabolites
    and cofactors).
    """

    node: MapNode
    score: float
    via: tuple[str, ...]
    accepted: bool
    vicinity_of: tuple[str, str] | None = None


@dataclass
class ReactionMapMatch:
    """Where a model reaction is on the map.

    `status` is ``on map`` if at least one candidate node is accepted,
    ``related`` if nodes with a common reaction ID or EC number were found but
    none matches the metabolites well enough, ``not on map`` otherwise and
    ``boundary`` for skipped boundary reactions.
    """

    reaction_id: str
    name: str
    status: str
    candidates: list[NodeMatch] = field(default_factory=list)

    @property
    def nodes(self) -> list[NodeMatch]:
        """The accepted candidates, best first."""
        return [c for c in self.candidates if c.accepted]


class MapMatcher:
    """Finds the nodes of a :class:`BrendaMap` that show the reactions of a model.

    Candidate nodes of a model reaction are found by reaction IDs (BRENDA
    reaction IDs of the node, aligned in BKMS-react with KEGG/MetaCyc
    reactions of the model reaction, including cross-references found in
    BiGGr), by EC numbers (all EC numbers found for the reaction, made current;
    incomplete EC numbers of a node include the complete ones) and by common
    main metabolites. Each candidate is scored by a weighted Jaccard index of
    the metabolites (as for the selection of EC numbers, cofactors count
    half). Metabolites are compared by KEGG compound IDs and by their names
    (see :func:`name_key`), map metabolites get KEGG IDs from KEGG's and BiGGr's
    compound names.

    A candidate found by reaction ID or EC number is accepted if its score is
    at least `min_score`, a candidate found only by metabolites if its score
    is at least `min_score_unlinked`; in both cases only if its score is at
    most `max_gap` lower than the best score of the reaction. So a reaction
    drawn in several subsystems is found on all these nodes, but not on nodes
    of variants of it with other cofactors (e.g. ADP- instead of
    ATP-dependent phosphofructokinase).

    Optionally, the nodes next to found nodes are searched with lower
    thresholds (`vicinity`, see :meth:`search_vicinity`).

    Parameters
    ----------
    brenda_map : BrendaMap
    mapper : ECMapper, optional
        For the EC numbers and identifiers of the reactions; default: ``ECMapper()``.
    min_score, min_score_unlinked, max_gap : float
    vicinity : str, optional
        ``neighbours`` or ``propagate`` (:data:`VICINITY_MODES`) to search the
        vicinity of the found nodes in :meth:`match_model`; default: no search.
    min_score_vicinity, min_score_vicinity_unlinked : float
        Thresholds for nodes in the vicinity, with and without common reaction ID or EC number.
    cache_dir, refresh :
        For the resources that are created (default: like the mapper's sources).
    """

    def __init__(
        self,
        brenda_map: BrendaMap,
        mapper: ECMapper | None = None,
        min_score: float = 0.6,
        min_score_unlinked: float = 0.8,
        max_gap: float = 0.1,
        vicinity: str | None = None,
        min_score_vicinity: float = 0.4,
        min_score_vicinity_unlinked: float = 0.6,
        cache_dir=None,
        refresh: bool = False,
    ):
        if vicinity is not None and vicinity not in VICINITY_MODES:
            raise ValueError(f"unknown vicinity mode {vicinity!r}, use one of {', '.join(VICINITY_MODES)}")
        self.map = brenda_map
        self.mapper = mapper or ECMapper(cache_dir=cache_dir, refresh=refresh)
        if cache_dir is None and self.mapper.sources:
            cache_dir = next(iter(self.mapper.sources.values())).cache_dir
        sources = list(self.mapper.sources.values())
        self.bkms = next((s for s in sources if isinstance(s, BKMSSource)), None) or BKMSSource(
            cache_dir=cache_dir, refresh=refresh
        )
        comparison = self.mapper.disambiguator.comparison
        if comparison is None:
            kegg = next((s for s in sources if isinstance(s, KEGGSource)), None) or KEGGSource(cache_dir=cache_dir)
            comparison = MetaboliteComparison(
                kegg,
                KEGGReactionCompounds(cache_dir=cache_dir, refresh=refresh),
                BiGGrMetabolites(cache_dir=cache_dir, refresh=refresh),
                self.mapper.disambiguator.enzyme,
            )
        self.comparison = comparison
        self.names = KEGGCompoundNames(cache_dir=cache_dir, refresh=refresh)
        self.min_score = min_score
        self.min_score_unlinked = min_score_unlinked
        self.max_gap = max_gap
        self.vicinity = vicinity
        self.min_score_vicinity = min_score_vicinity
        self.min_score_vicinity_unlinked = min_score_vicinity_unlinked
        self._index = None

    @property
    def resources(self) -> list[CachedResource]:
        """All resources with cached files used for the matching."""
        resources = list(self.mapper.resources)
        for resource in [self.bkms, *self.comparison.resources, self.names]:
            if all(resource is not r for r in resources):
                resources.append(resource)
        return resources

    def reset(self) -> None:
        """Forget the indexes (after an update of the resources)."""
        self._index = None

    # --- metabolites ------------------------------------------------------------------------

    def _name_index(self) -> dict[str, set[str]]:
        index: dict[str, set[str]] = {}
        for compound, names in self.names.data.items():
            for name in names:
                index.setdefault(name_key(name), set()).add(compound)
        bigg = self.comparison.bigg_metabolites
        for xrefs in (bigg.data.values() if bigg is not None else ()):
            for name in xrefs.get("name", ()):
                index.setdefault(name_key(name), set()).update(xrefs.get("kegg.compound", ()))
        return index

    def ligand_keys(self, names: Iterable[str]) -> frozenset[str]:
        """Keys (KEGG compound IDs and ``name:`` keys) of a map metabolite; empty for protons."""
        keys = {name_key(n) for n in expand_names(names)} - _PROTON_NAMES - {""}
        if not keys:
            return frozenset()
        compounds = set().union(*(self._index["names"].get(k, set()) for k in keys))
        return _keys(compounds - {PROTON}, keys)

    def metabolite_keys(self, metabolite) -> frozenset[str]:
        """Keys (KEGG compound IDs and ``name:`` keys) of a cobrapy metabolite; empty for protons."""
        compounds = self.comparison.metabolite_compounds(metabolite)
        names = {metabolite.name} if metabolite.name else set()
        bigg = self.comparison.bigg_metabolites
        if bigg is not None:
            for bigg_id in self.comparison.bigg_ids(metabolite):
                names |= bigg.xrefs(bigg_id).get("name", set())
        for compound in compounds:
            names.update(self.names.names(compound))
        keys = {name_key(n) for n in names}
        if not compounds and (keys & _PROTON_NAMES or PROTON in annotation.annotation_values(metabolite.annotation, "kegg.compound")):
            return frozenset()
        return _keys(compounds, keys)

    def reaction_keys(self, reaction) -> list[frozenset[str]]:
        """Keys of the metabolites of a cobrapy reaction, one set per metabolite, without duplicates."""
        result = []
        for metabolite in reaction.metabolites:
            keys = self.metabolite_keys(metabolite)
            if keys and keys not in result:
                result.append(keys)
        return result

    def is_transport(self, reaction) -> bool:
        """Whether a metabolite (other than protons) is on both sides of the reaction, e.g. in two compartments."""
        sides = ({}, {})
        for metabolite, coefficient in reaction.metabolites.items():
            keys = self.metabolite_keys(metabolite)
            if keys:
                sides[coefficient > 0][keys] = True
        return bool(sides[0].keys() & sides[1].keys())

    @staticmethod
    def score(left: list[frozenset[str]], right: list[frozenset[str]], currency_weight: float = 0.5) -> float:
        """Weighted Jaccard index of two lists of metabolites (sets of keys); currency metabolites count less."""
        matched_right, weight_matched, weight_union = set(), 0.0, 0.0
        for keys in left:
            weight = currency_weight if keys & CURRENCY_KEYS else 1.0
            weight_union += weight
            hits = [i for i, other in enumerate(right) if keys & other]
            if hits:
                weight_matched += weight
                matched_right.update(hits)
        for i, other in enumerate(right):
            if i not in matched_right:
                weight_union += currency_weight if other & CURRENCY_KEYS else 1.0
        return weight_matched / weight_union if weight_union else 0.0

    # --- index ------------------------------------------------------------------------------

    def _build(self) -> None:
        self._index = {"names": self._name_index()}
        enzyme = self.mapper.disambiguator.enzyme
        aligned_kegg = {
            node.node_id: set().union(
                set(), *(self.bkms.brenda_reaction(br).get("kegg.reaction", set()) for br in node.brenda_reaction_ids)
            )
            for node in self.map.nodes.values()
        }
        inferred = self._infer_compounds(aligned_kegg)
        by_ec, by_partial_ec, by_xref, by_main = {}, {}, {}, {}
        node_keys, node_ecs = {}, {}
        for node in self.map.nodes.values():
            node_ecs[node.node_id] = set()
            for raw in node.ec_numbers:
                ec = normalize_ec(raw)
                if ec is None:
                    continue
                for current in (enzyme.current(ec) if enzyme is not None else {ec}) or {ec}:
                    target = by_ec if is_complete(current) else by_partial_ec
                    target.setdefault(current, set()).add(node.node_id)
                    node_ecs[node.node_id].add(current)
            for br in node.brenda_reaction_ids:
                for namespace, ids in self.bkms.brenda_reaction(br).items():
                    for identifier in ids:
                        by_xref.setdefault((namespace, identifier), set()).add(node.node_id)
            keys = []
            for metabolite in node.metabolites:
                k = self.ligand_keys(metabolite.names)
                if k and metabolite.ligand_id in inferred:
                    k = k | {inferred[metabolite.ligand_id]}
                if k and k not in keys:
                    keys.append(k)
                if k and metabolite.role == "main" and not k & CURRENCY_KEYS:
                    for key in k:
                        by_main.setdefault(key, set()).add(node.node_id)
            node_keys[node.node_id] = keys
        self._index.update(
            by_ec=by_ec,
            by_partial_ec=by_partial_ec,
            by_xref=by_xref,
            by_main=by_main,
            keys=node_keys,
            ecs=node_ecs,
            inferred=inferred,
            places={},
        )

    def _infer_compounds(self, aligned_kegg: dict[str, set[str]]) -> dict[int, str]:
        """KEGG compounds of map metabolites whose names are not known, {ligand ID: compound}.

        If the BRENDA reaction of a node is aligned with a KEGG reaction and
        exactly one metabolite of the node has no KEGG compound and exactly one
        compound of the KEGG reaction is not matched by the other metabolites,
        they are taken to be the same. The compound found most often for a
        ligand ID is used for all nodes with this ligand.
        """
        votes: dict[int, Counter] = {}
        for node in self.map.nodes.values():
            if not aligned_kegg[node.node_id]:
                continue
            keyed = [(m, self.ligand_keys(m.names)) for m in node.metabolites]
            keyed = [(m, k) for m, k in keyed if k]
            unresolved = [m for m, k in keyed if not any(_KEGG_ID.fullmatch(x) for x in k)]
            if len(unresolved) != 1 or unresolved[0].ligand_id is None:
                continue
            known = set().union(set(), *(k for _, k in keyed))
            for reaction in aligned_kegg[node.node_id]:
                free = set(self.comparison.compounds.compounds(reaction)) - {PROTON} - known
                if len(free) == 1:
                    votes.setdefault(unresolved[0].ligand_id, Counter())[free.pop()] += 1
        return {ligand: counter.most_common(1)[0][0] for ligand, counter in votes.items()}

    def node_keys(self, node_id: str) -> list[frozenset[str]]:
        if self._index is None:
            self._build()
        return self._index["keys"][node_id]

    def place_keys(self, place: str) -> frozenset[str]:
        """Keys of the main metabolite drawn at a place of the map (see :meth:`BrendaMap.neighbours`)."""
        if self._index is None:
            self._build()
        cache = self._index["places"]
        if place not in cache:
            metabolites = self.map.place_metabolites(place)
            keys = self.ligand_keys([n for m in metabolites for n in m.names])
            inferred = {self._index["inferred"].get(m.ligand_id) for m in metabolites} - {None}
            cache[place] = keys | inferred if keys else keys
        return cache[place]

    def compatible_ecs(self, reaction_ecs: Iterable[str], node_id: str) -> bool:
        """Whether the EC numbers of a reaction do not contradict those of a node.

        They contradict if the reaction has complete EC numbers, the node EC
        numbers, and no pair of them has the same first three levels (``-`` in
        an EC number of the node agrees with all).
        """
        if self._index is None:
            self._build()
        complete = [ec.split(".")[:3] for ec in reaction_ecs if is_complete(ec)]
        node = [ec.split(".")[:3] for ec in self._index["ecs"][node_id]]
        if not complete or not node:
            return True
        return any(all(b == "-" or a == b for a, b in zip(r, n)) for r in complete for n in node)

    # --- matching ---------------------------------------------------------------------------

    def reaction_ecs(self, reaction) -> set[str]:
        """All EC numbers of a reaction: from the model and from all resources, made current."""
        ecs = set(annotation.annotated_ecs(reaction)) | set(annotation.notes_ecs(reaction))
        ecs |= {e.ec for e in self.mapper.query(annotation.xrefs(reaction))}
        return set(self.mapper.disambiguator.current(ecs))

    def match_reaction(self, reaction, boundary: bool = False) -> ReactionMapMatch:
        """Find the map nodes of a cobrapy reaction."""
        if self._index is None:
            self._build()
        index = self._index
        if boundary:
            return ReactionMapMatch(reaction.id, reaction.name, "boundary")
        via: dict[str, set[str]] = {}
        identifiers = self.mapper.expand_xrefs(annotation.xrefs(reaction))
        for namespace in ("kegg.reaction", "metacyc.reaction"):
            for identifier in identifiers.get(namespace, {}):
                for node_id in index["by_xref"].get((namespace, identifier), ()):
                    via.setdefault(node_id, set()).add("reaction-id")
        for ec in self.reaction_ecs(reaction):
            nodes = set(index["by_ec"].get(ec, ()))
            for partial, partial_nodes in index["by_partial_ec"].items():
                if covers(partial, ec):
                    nodes |= partial_nodes
            for node_id in nodes:
                via.setdefault(node_id, set()).add("ec")
        keys = self.reaction_keys(reaction)
        if not self.is_transport(reaction):
            # transporters are not drawn, their metabolites would match enzymes that convert them
            shared = Counter(
                node_id
                for k in keys
                if not k & CURRENCY_KEYS
                for node_id in set().union(*(index["by_main"].get(key, set()) for key in k))
            )
            for node_id, count in shared.items():
                if count >= 2:
                    via.setdefault(node_id, set()).add("metabolites")

        candidates = []
        for node_id, how in via.items():
            score = self.score(keys, index["keys"][node_id])
            linked = bool(how & {"reaction-id", "ec"})
            accepted = score >= (self.min_score if linked else self.min_score_unlinked)
            order = tuple(v for v in ("reaction-id", "ec", "metabolites") if v in how)
            candidates.append(NodeMatch(self.map.nodes[node_id], score, order, accepted))
        best = max((c.score for c in candidates if c.accepted), default=None)
        for c in candidates:
            if c.accepted and c.score < best - self.max_gap:
                c.accepted = False
        candidates.sort(key=_candidate_order)
        if any(c.accepted for c in candidates):
            status = "on map"
        elif any(set(c.via) & {"reaction-id", "ec"} for c in candidates):
            status = "related"
        else:
            status = "not on map"
        return ReactionMapMatch(reaction.id, reaction.name, status, candidates)

    def match_model(self, model, reactions: Iterable | None = None, skip_boundary: bool = True) -> dict[str, ReactionMapMatch]:
        """Find the map nodes of the reactions of a cobrapy model (default: all) as {reaction ID: match}.

        If :attr:`vicinity` is set, the vicinity of the found nodes is searched afterwards.
        """
        boundary = {r.id for r in model.boundary} if skip_boundary else set()
        reactions = [
            model.reactions.get_by_id(r) if isinstance(r, str) else r
            for r in (model.reactions if reactions is None else reactions)
        ]
        result = {r.id: self.match_reaction(r, r.id in boundary) for r in reactions}
        if self.vicinity:
            self.search_vicinity(reactions, result, self.vicinity)
        return result

    # --- vicinity ---------------------------------------------------------------------------

    def search_vicinity(
        self, reactions: Iterable, matches: dict[str, ReactionMapMatch], mode: str = "neighbours"
    ) -> list[tuple[str, str]]:
        """Accept further nodes for model reactions next to the nodes already found.

        Two enzyme nodes of the map are neighbours if they share a drawn main
        metabolite (or metabolite nodes connected by a link) that is not a
        currency metabolite. For a node N' next to node N of model reaction R,
        the candidates are the model reactions R' that have this metabolite,
        too, if R also has it: the connection on the map is then also one in the
        model. Transporters are no candidates. A candidate (R', N') is accepted
        if

        - its score is at least `min_score_vicinity` if R' and N' have a common
          reaction ID or EC number, otherwise at least
          `min_score_vicinity_unlinked`, and their EC numbers do not
          contradict each other (:meth:`compatible_ecs`), and
        - its score is at most `max_gap` lower than that of the best candidate
          or model reaction of N', and than that of the best node of R'.

        Modes:

        ``neighbours``
            Starting from all found nodes, the neighbours without model
            reactions are searched once.
        ``propagate``
            Starting from the unambiguous matches (reactions with only one
            accepted node, which has a common reaction ID or EC number), all
            neighbours are searched. The accepted nodes and the found nodes
            that were reached are the starting points of the next round, until
            nothing new is found. Found nodes that are not reached stay.

        `matches` (from :meth:`match_reaction` for `reactions`) are changed in
        place: accepted nodes get ``vicinity`` in `via` and the match they
        are next to as `vicinity_of`, reactions with accepted nodes are ``on
        map``. Returns the new matches as (reaction ID, node ID).
        """
        if mode not in VICINITY_MODES:
            raise ValueError(f"unknown vicinity mode {mode!r}, use one of {', '.join(VICINITY_MODES)}")
        if self._index is None:
            self._build()
        reactions = {r.id: r for r in reactions if r.id in matches and matches[r.id].status != "boundary"}
        keys = {rid: self.reaction_keys(r) for rid, r in reactions.items()}
        by_key: dict[str, set[str]] = {}
        for rid, reaction in reactions.items():
            if self.is_transport(reaction):
                continue
            for k in keys[rid]:
                if not k & CURRENCY_KEYS:
                    for key in k:
                        by_key.setdefault(key, set()).add(rid)

        accepted = {(rid, c.node.node_id) for rid in reactions for c in matches[rid].nodes}
        if mode == "propagate":
            n_nodes = Counter(rid for rid, _ in accepted)
            anchors = {
                (rid, c.node.node_id)
                for rid in reactions
                for c in matches[rid].nodes
                if n_nodes[rid] == 1 and set(c.via) & {"reaction-id", "ec"}
            }
        else:
            anchors = set(accepted)
        identified = {node_id for _, node_id in accepted}
        ecs: dict[str, set[str]] = {}
        frontier, new = set(anchors), []
        while frontier:
            candidates: dict[tuple[str, str], tuple[str, str]] = {}
            for rid, node_id in sorted(frontier):
                for other, place in self.map.neighbours(node_id):
                    if mode == "neighbours" and other in identified:
                        continue
                    shared = self.place_keys(place)
                    if not shared or shared & CURRENCY_KEYS or not any(shared & k for k in keys[rid]):
                        continue
                    for rid2 in sorted(set().union(*(by_key.get(key, set()) for key in shared))):
                        candidates.setdefault((rid2, other), (rid, node_id))
            reached = {pair for pair in candidates if pair in accepted}
            chosen = self._select_vicinity([p for p in candidates if p not in accepted], matches, reactions, keys, ecs)
            for (rid, node_id), score in chosen.items():
                match = matches[rid]
                candidate = next((c for c in match.candidates if c.node.node_id == node_id), None)
                if candidate is None:
                    candidate = NodeMatch(self.map.nodes[node_id], score, (), False)
                    match.candidates.append(candidate)
                candidate.accepted = True
                candidate.via += ("vicinity",)
                candidate.vicinity_of = candidates[(rid, node_id)]
                match.status = "on map"
                match.candidates.sort(key=_candidate_order)
                new.append((rid, node_id))
            accepted |= chosen.keys()
            if mode == "neighbours":
                break
            frontier = (reached | chosen.keys()) - anchors
            anchors |= frontier
        return new

    def _select_vicinity(self, pairs, matches, reactions, keys, ecs) -> dict[tuple[str, str], float]:
        """The accepted candidates (reaction ID, node ID) of the vicinity with their scores."""
        eligible = {}
        for rid, node_id in pairs:
            existing = next((c for c in matches[rid].candidates if c.node.node_id == node_id), None)
            if existing is not None:
                score, linked = existing.score, bool(set(existing.via) & {"reaction-id", "ec"})
            else:
                score, linked = self.score(keys[rid], self.node_keys(node_id)), False
            if linked:
                ok = score >= self.min_score_vicinity
            else:
                if rid not in ecs:
                    ecs[rid] = self.reaction_ecs(reactions[rid])
                ok = score >= self.min_score_vicinity_unlinked and self.compatible_ecs(ecs[rid], node_id)
            if ok:
                eligible[(rid, node_id)] = score
        best_node: dict[str, float] = {}
        best_reaction: dict[str, float] = {}
        accepted = [(rid, c.node.node_id, c.score) for rid in reactions for c in matches[rid].nodes]
        for rid, node_id, score in accepted + [(r, n, s) for (r, n), s in eligible.items()]:
            best_node[node_id] = max(best_node.get(node_id, 0.0), score)
            best_reaction[rid] = max(best_reaction.get(rid, 0.0), score)
        return {
            (rid, node_id): score
            for (rid, node_id), score in eligible.items()
            if score >= best_node[node_id] - self.max_gap and score >= best_reaction[rid] - self.max_gap
        }


def _candidate_order(c: NodeMatch):
    return (-c.accepted, -c.score, -len(c.via), c.node.node_id)


# --- output ---------------------------------------------------------------------------------


def _node_columns(node: MapNode | None) -> dict:
    return {
        "node_id": node.node_id if node else "",
        "subsystem": node.subsystem if node else "",
        "ec_label": node.ec_label if node else "",
        "enzyme_name": node.enzyme_name if node else "",
        "x": node.x if node else None,
        "y": node.y if node else None,
    }


def map_matches_to_dataframe(matches: dict[str, ReactionMapMatch], all_candidates: bool = False) -> pd.DataFrame:
    """One row per model reaction and node.

    Reactions on the map get one row per accepted node, related reactions one
    row for their best candidate node, others one row without node. With
    `all_candidates`, all candidate nodes are listed. For nodes found in the
    vicinity of another match, `vicinity_of` is that match ("reaction@node").
    """
    rows = []
    for m in matches.values():
        if all_candidates:
            shown = m.candidates
        elif m.status == "on map":
            shown = m.nodes
        else:
            shown = [c for c in m.candidates if set(c.via) & {"reaction-id", "ec"}][:1]
        for c in shown or [None]:
            rows.append(
                {
                    "reaction_id": m.reaction_id,
                    "reaction_name": m.name,
                    "status": m.status,
                    **_node_columns(c.node if c else None),
                    "score": round(c.score, 3) if c else None,
                    "via": ";".join(c.via) if c else "",
                    "vicinity_of": "@".join(c.vicinity_of) if c and c.vicinity_of else "",
                    "accepted": c.accepted if c else False,
                }
            )
    columns = ["reaction_id", "reaction_name", "status", "node_id", "subsystem", "ec_label", "enzyme_name",
               "score", "via", "vicinity_of", "accepted", "x", "y"]
    return pd.DataFrame(rows, columns=columns)


def node_coverage(brenda_map: BrendaMap, matches: dict[str, ReactionMapMatch]) -> pd.DataFrame:
    """One row per enzyme node of the map with the model reactions shown by it."""
    reactions: dict[str, list[str]] = {}
    for m in matches.values():
        for c in m.nodes:
            reactions.setdefault(c.node.node_id, []).append(m.reaction_id)
    rows = [
        {**_node_columns(node), "model_reactions": ";".join(reactions.get(node.node_id, []))}
        for node in brenda_map.nodes.values()
    ]
    columns = ["node_id", "subsystem", "ec_label", "enzyme_name", "x", "y", "model_reactions"]
    return pd.DataFrame(rows, columns=columns)


def _only_vicinity(matches: dict[str, ReactionMapMatch]) -> tuple[set[str], set[str]]:
    """Reactions and nodes that are on the map only by matches in the vicinity of others."""
    regular_nodes = {c.node.node_id for m in matches.values() for c in m.nodes if "vicinity" not in c.via}
    reactions = {m.reaction_id for m in matches.values() if m.nodes and all("vicinity" in c.via for c in m.nodes)}
    nodes = {c.node.node_id for m in matches.values() for c in m.nodes} - regular_nodes
    return reactions, nodes


def map_summary(brenda_map: BrendaMap, matches: dict[str, ReactionMapMatch]) -> dict:
    """Reactions per status, and covered nodes and subsystems of the map.

    `vicinity_reactions` and `vicinity_nodes` count the reactions and nodes
    that are on the map only by matches found in the vicinity of others.
    """
    statuses = Counter(m.status for m in matches.values())
    nodes = {c.node.node_id for m in matches.values() for c in m.nodes}
    subsystems = Counter(brenda_map.nodes[n].subsystem for n in nodes)
    vicinity_reactions, vicinity_nodes = _only_vicinity(matches)
    return {
        "reactions": len(matches),
        "status": {s: statuses.get(s, 0) for s in MAP_STATUSES},
        "nodes": len(brenda_map),
        "nodes_with_model_reactions": len(nodes),
        "subsystems_with_model_reactions": len(subsystems),
        "nodes_per_subsystem": dict(subsystems.most_common()),
        "vicinity_reactions": len(vicinity_reactions),
        "vicinity_nodes": len(vicinity_nodes),
    }


_KEGG_ID = re.compile(r"C\d{5}")

_SVG_STYLE = """
g.nodeEN.ecm-on-map rect{fill:rgb(110, 200, 100);stroke:rgb(50, 130, 45);stroke-width:3;}
g.nodeEN.ecm-vicinity rect{fill:rgb(80, 200, 200);stroke:rgb(30, 120, 130);stroke-width:3;}
g.nodeEN.ecm-related rect{fill:rgb(255, 225, 110);stroke:rgb(190, 160, 40);}
"""

_SVG_NODE = re.compile(r'<g id="(node\d+)" class="([^"]*\bnodeEN\b[^"]*)"([^>]*)>(.*?)</g>', re.DOTALL)


def highlight_svg(
    svg: str | os.PathLike, output: str | os.PathLike, matches: dict[str, ReactionMapMatch], related: bool = True
) -> int:
    """Write a copy of the map SVG in which the nodes of model reactions are highlighted.

    Nodes with accepted model reactions get the class ``ecm-on-map`` (green),
    or ``ecm-vicinity`` (cyan) if all were found in the vicinity of other
    matches, with `related` the best candidate nodes of related reactions
    ``ecm-related`` (yellow). The IDs of the model reactions are added to the
    tooltip (``title``) of the nodes. Returns the number of highlighted nodes.
    """
    on_map: dict[str, list[str]] = {}
    near: dict[str, list[str]] = {}
    for m in matches.values():
        for c in m.nodes:
            on_map.setdefault(c.node.node_id, []).append(m.reaction_id)
        if related and m.status == "related":
            best = [c for c in m.candidates if set(c.via) & {"reaction-id", "ec"}][:1]
            for c in best:
                near.setdefault(c.node.node_id, []).append(m.reaction_id)
    vicinity = _only_vicinity(matches)[1]
    text = Path(svg).read_text(encoding="utf-8")
    highlighted = 0

    def replace(match: re.Match) -> str:
        nonlocal highlighted
        node_id, classes, attributes, content = match.groups()
        if node_id in vicinity:
            css, label, ids = "ecm-vicinity", "model (vicinity)", on_map[node_id]
        elif node_id in on_map:
            css, label, ids = "ecm-on-map", "model", on_map[node_id]
        elif node_id in near:
            css, label, ids = "ecm-related", "related model reactions", near[node_id]
        else:
            return match.group(0)
        highlighted += 1
        note = html.escape(f"{label}: {', '.join(ids)}")
        if "<title>" in content:
            content = content.replace("</title>", f" | {note}</title>", 1)
        else:
            content += f"<title>{note}</title>"
        return f'<g id="{node_id}" class="{classes} {css}"{attributes}>{content}</g>'

    text = _SVG_NODE.sub(replace, text)
    text = text.replace("</style>", _SVG_STYLE + "</style>", 1)
    Path(output).write_text(text, encoding="utf-8")
    return highlighted
