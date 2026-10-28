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
    """The enzyme nodes of a BRENDA pathway map."""

    def __init__(self, nodes: Iterable[MapNode]):
        self.nodes = {n.node_id: n for n in nodes}

    def __len__(self) -> int:
        return len(self.nodes)

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
        return cls(nodes)


@dataclass
class NodeMatch:
    """A map node that is a candidate for a model reaction.

    `via` lists how the node was found: ``reaction-id`` (a BRENDA reaction of
    the node is aligned in BKMS-react with a KEGG or MetaCyc reaction of the
    model reaction), ``ec`` (common EC number) and ``metabolites`` (at least two
    common main metabolites). `score` compares the metabolites of the model
    reaction with those of the node (main metabolites and cofactors).
    """

    node: MapNode
    score: float
    via: tuple[str, ...]
    accepted: bool


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

    Parameters
    ----------
    brenda_map : BrendaMap
    mapper : ECMapper, optional
        For the EC numbers and identifiers of the reactions; default: ``ECMapper()``.
    min_score, min_score_unlinked, max_gap : float
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
        cache_dir=None,
        refresh: bool = False,
    ):
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
        node_keys = {}
        for node in self.map.nodes.values():
            for raw in node.ec_numbers:
                ec = normalize_ec(raw)
                if ec is None:
                    continue
                for current in (enzyme.current(ec) if enzyme is not None else {ec}) or {ec}:
                    target = by_ec if is_complete(current) else by_partial_ec
                    target.setdefault(current, set()).add(node.node_id)
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
            by_ec=by_ec, by_partial_ec=by_partial_ec, by_xref=by_xref, by_main=by_main, keys=node_keys, inferred=inferred
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
        candidates.sort(key=lambda c: (-c.accepted, -c.score, -len(c.via), c.node.node_id))
        if any(c.accepted for c in candidates):
            status = "on map"
        elif any(set(c.via) & {"reaction-id", "ec"} for c in candidates):
            status = "related"
        else:
            status = "not on map"
        return ReactionMapMatch(reaction.id, reaction.name, status, candidates)

    def match_model(self, model, reactions: Iterable | None = None, skip_boundary: bool = True) -> dict[str, ReactionMapMatch]:
        """Find the map nodes of the reactions of a cobrapy model (default: all) as {reaction ID: match}."""
        boundary = {r.id for r in model.boundary} if skip_boundary else set()
        reactions = model.reactions if reactions is None else reactions
        result = {}
        for reaction in reactions:
            if isinstance(reaction, str):
                reaction = model.reactions.get_by_id(reaction)
            result[reaction.id] = self.match_reaction(reaction, reaction.id in boundary)
        return result


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
    `all_candidates`, all candidate nodes are listed.
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
                    "accepted": c.accepted if c else False,
                }
            )
    columns = ["reaction_id", "reaction_name", "status", "node_id", "subsystem", "ec_label", "enzyme_name",
               "score", "via", "accepted", "x", "y"]
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


def map_summary(brenda_map: BrendaMap, matches: dict[str, ReactionMapMatch]) -> dict:
    """Reactions per status, and covered nodes and subsystems of the map."""
    statuses = Counter(m.status for m in matches.values())
    nodes = {c.node.node_id for m in matches.values() for c in m.nodes}
    subsystems = Counter(brenda_map.nodes[n].subsystem for n in nodes)
    return {
        "reactions": len(matches),
        "status": {s: statuses.get(s, 0) for s in MAP_STATUSES},
        "nodes": len(brenda_map),
        "nodes_with_model_reactions": len(nodes),
        "subsystems_with_model_reactions": len(subsystems),
        "nodes_per_subsystem": dict(subsystems.most_common()),
    }


_KEGG_ID = re.compile(r"C\d{5}")

_SVG_STYLE = """
g.nodeEN.ecm-on-map rect{fill:rgb(110, 200, 100);stroke:rgb(50, 130, 45);stroke-width:3;}
g.nodeEN.ecm-related rect{fill:rgb(255, 225, 110);stroke:rgb(190, 160, 40);}
"""

_SVG_NODE = re.compile(r'<g id="(node\d+)" class="([^"]*\bnodeEN\b[^"]*)"([^>]*)>(.*?)</g>', re.DOTALL)


def highlight_svg(
    svg: str | os.PathLike, output: str | os.PathLike, matches: dict[str, ReactionMapMatch], related: bool = True
) -> int:
    """Write a copy of the map SVG in which the nodes of model reactions are highlighted.

    Nodes with accepted model reactions get the class ``ecm-on-map`` (green),
    with `related` the best candidate nodes of related reactions
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
    text = Path(svg).read_text(encoding="utf-8")
    highlighted = 0

    def replace(match: re.Match) -> str:
        nonlocal highlighted
        node_id, classes, attributes, content = match.groups()
        if node_id in on_map:
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
