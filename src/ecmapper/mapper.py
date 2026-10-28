"""Mapping of model reactions to EC numbers."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable

import pandas as pd

from . import annotation
from .ec import is_complete
from .sources import ECSource, default_sources

#: values of :attr:`ReactionECMapping.origin`
ORIGINS = ("annotation", "notes", "resources", "none", "boundary")


@dataclass(frozen=True)
class ECEvidence:
    """One EC number for a reaction and where it came from.

    Attributes
    ----------
    ec : str
        The EC number.
    source : str
        ``annotation``, ``notes`` or the name of the resource that was queried.
    identifier : str
        The identifier that was looked up in `source` (empty for annotation and notes).
    via : str
        For identifiers that were not in the model but were found as
        cross-reference in another resource: that resource's identifier as
        ``namespace:identifier``, e.g. ``bigg.reaction:PGI``.
    """

    ec: str
    source: str
    identifier: str = ""
    via: str = ""


@dataclass
class ReactionECMapping:
    """The EC numbers found for one reaction.

    Attributes
    ----------
    reaction_id : str
    origin : str
        ``annotation`` or ``notes`` if the model already had EC numbers,
        ``resources`` if they were found by querying the resources, ``none`` if
        nothing was found and ``boundary`` for skipped boundary reactions.
    evidence : list of ECEvidence
    xrefs : dict
        The cross-references of the reaction in the model, {namespace: [identifiers]}.
    """

    reaction_id: str
    origin: str
    evidence: list[ECEvidence] = field(default_factory=list)
    xrefs: dict[str, list[str]] = field(default_factory=dict)

    @property
    def ec_numbers(self) -> list[str]:
        """All EC numbers, sorted."""
        return sorted({e.ec for e in self.evidence})

    @property
    def ambiguous(self) -> bool:
        """True if more than one EC number was found."""
        return len(self.ec_numbers) > 1

    @property
    def sources(self) -> list[str]:
        """The sources that contributed EC numbers, in order of their first contribution."""
        return list(dict.fromkeys(e.source for e in self.evidence))

    def sources_of(self, ec: str) -> list[str]:
        """The sources that support `ec`."""
        return list(dict.fromkeys(e.source for e in self.evidence if e.ec == ec))


class ECMapper:
    """Maps reactions of a cobrapy model to EC numbers.

    For every reaction the EC numbers in the annotation are used. If there are
    none, EC numbers in COBRA style notes (``EC Number: ...``) are used. If
    there are none either, the cross-references in the annotation (BiGG, KEGG,
    ModelSEED, Rhea, MetaCyc/BioCyc) are looked up in the resources and the
    union of the EC numbers found is returned, together with their evidence.

    Parameters
    ----------
    sources : iterable of ECSource or str, optional
        The resources to query, instances or names (see :data:`ecmapper.sources.SOURCES`).
        Default: all resources.
    use_notes : bool
        Fall back to EC numbers in the notes of a reaction.
    always_query : bool
        Query the resources also for reactions that already have EC numbers.
    skip_boundary : bool
        Do not map boundary (exchange, demand, sink) reactions.
    follow_xrefs : bool
        Also look up identifiers that a resource lists as cross-references of
        the reaction (one step, e.g. BiGG ID -> KEGG reaction in BiGGr -> EC).
    cache_dir, refresh :
        Passed to the default sources, see :class:`ecmapper.sources.ECSource`.
    """

    def __init__(
        self,
        sources: Iterable[ECSource | str] | None = None,
        *,
        use_notes: bool = True,
        always_query: bool = False,
        skip_boundary: bool = True,
        follow_xrefs: bool = True,
        cache_dir=None,
        refresh: bool = False,
    ):
        if sources is None:
            sources = default_sources(cache_dir=cache_dir, refresh=refresh)
        else:
            sources = list(sources)
            names = [s for s in sources if isinstance(s, str)]
            created = iter(default_sources(names, cache_dir=cache_dir, refresh=refresh))
            sources = [next(created) if isinstance(s, str) else s for s in sources]
        self.sources: dict[str, ECSource] = {}
        for source in sources:
            if source.namespace in self.sources:
                raise ValueError(f"more than one source for namespace {source.namespace}")
            self.sources[source.namespace] = source
        self.use_notes = use_notes
        self.always_query = always_query
        self.skip_boundary = skip_boundary
        self.follow_xrefs = follow_xrefs

    def map_reaction(self, reaction, boundary: bool = False) -> ReactionECMapping:
        """Map one cobrapy reaction, `boundary` marks it as boundary reaction."""
        xrefs = annotation.xrefs(reaction)
        if boundary and self.skip_boundary:
            return ReactionECMapping(reaction.id, "boundary", xrefs=xrefs)
        origin = "none"
        evidence = [ECEvidence(ec, "annotation") for ec in annotation.annotated_ecs(reaction)]
        if evidence:
            origin = "annotation"
        elif self.use_notes:
            evidence = [ECEvidence(ec, "notes") for ec in annotation.notes_ecs(reaction)]
            if evidence:
                origin = "notes"
        if not evidence or self.always_query:
            found = self.query(xrefs)
            evidence += found
            if found and origin == "none":
                origin = "resources"
        return ReactionECMapping(reaction.id, origin, evidence, xrefs)

    def map_model(self, model, reactions: Iterable | None = None) -> dict[str, ReactionECMapping]:
        """Map the reactions of a cobrapy model (default: all) as {reaction ID: mapping}."""
        boundary = {r.id for r in model.boundary}
        reactions = model.reactions if reactions is None else reactions
        result = {}
        for reaction in reactions:
            if isinstance(reaction, str):
                reaction = model.reactions.get_by_id(reaction)
            result[reaction.id] = self.map_reaction(reaction, reaction.id in boundary)
        return result

    def query(self, xrefs: dict[str, list[str]]) -> list[ECEvidence]:
        """Look up the cross-references {namespace: [identifiers]} in the resources."""
        evidence = []
        queried = set()
        for namespace, identifiers in xrefs.items():
            source = self.sources.get(namespace)
            if source is None:
                continue
            for identifier in identifiers:
                queried.add((namespace, identifier))
                evidence += [ECEvidence(ec, source.name, identifier) for ec in sorted(source.lookup(identifier))]
        if not self.follow_xrefs:
            return evidence
        for namespace, identifiers in xrefs.items():
            source = self.sources.get(namespace)
            if source is None:
                continue
            for identifier in identifiers:
                for xref_namespace, xref_ids in sorted(source.xrefs(identifier).items()):
                    xref_source = self.sources.get(xref_namespace)
                    if xref_source is None:
                        continue
                    for xref_id in sorted(xref_ids):
                        if (xref_namespace, xref_id) in queried:
                            continue
                        queried.add((xref_namespace, xref_id))
                        evidence += [
                            ECEvidence(ec, xref_source.name, xref_id, f"{namespace}:{identifier}")
                            for ec in sorted(xref_source.lookup(xref_id))
                        ]
        return evidence


def map_model(model, **kwargs) -> dict[str, ReactionECMapping]:
    """Shortcut for ``ECMapper(**kwargs).map_model(model)``."""
    return ECMapper(**kwargs).map_model(model)


def _join(values: Iterable[str]) -> str:
    return ";".join(values)


def to_dataframe(mappings: dict[str, ReactionECMapping] | Iterable[ReactionECMapping], long: bool = False) -> pd.DataFrame:
    """Tabulate mappings.

    The default (wide) table has one row per reaction with the EC numbers
    separated by ``;``. The long table has one row per reaction, EC number and
    evidence, reactions without EC number get one row with empty EC.
    """
    if isinstance(mappings, dict):
        mappings = mappings.values()
    rows = []
    for m in mappings:
        if long:
            for e in m.evidence or [ECEvidence("", "")]:
                rows.append(
                    {
                        "reaction_id": m.reaction_id,
                        "origin": m.origin,
                        "ec": e.ec,
                        "complete": is_complete(e.ec) if e.ec else None,
                        "source": e.source,
                        "identifier": e.identifier,
                        "via": e.via,
                    }
                )
        else:
            rows.append(
                {
                    "reaction_id": m.reaction_id,
                    "origin": m.origin,
                    "ec_numbers": _join(m.ec_numbers),
                    "n_ec": len(m.ec_numbers),
                    "sources": _join(m.sources),
                    "xrefs": _join(f"{ns}:{i}" for ns, ids in m.xrefs.items() for i in ids),
                }
            )
    columns = (
        ["reaction_id", "origin", "ec", "complete", "source", "identifier", "via"]
        if long
        else ["reaction_id", "origin", "ec_numbers", "n_ec", "sources", "xrefs"]
    )
    return pd.DataFrame(rows, columns=columns)


def summary(mappings: dict[str, ReactionECMapping]) -> dict:
    """Counts of reactions per origin, of ambiguous reactions and per source."""
    values = list(mappings.values())
    origins = Counter(m.origin for m in values)
    per_source = Counter(s for m in values for s in m.sources)
    return {
        "reactions": len(values),
        "origin": {o: origins.get(o, 0) for o in ORIGINS},
        "with_ec": sum(1 for m in values if m.evidence),
        "ambiguous": sum(1 for m in values if m.ambiguous),
        "reactions_per_source": dict(per_source.most_common()),
    }


def annotate_model(model, mappings: dict[str, ReactionECMapping], overwrite: bool = False) -> int:
    """Write EC numbers found in the resources to ``reaction.annotation["ec-code"]``.

    Only reactions with origin ``resources`` are changed, existing ``ec-code``
    entries are kept unless `overwrite`. Returns the number of changed reactions.
    """
    changed = 0
    for reaction_id, m in mappings.items():
        if m.origin != "resources":
            continue
        reaction = model.reactions.get_by_id(reaction_id)
        if "ec-code" in reaction.annotation and not overwrite:
            continue
        reaction.annotation["ec-code"] = m.ec_numbers
        changed += 1
    return changed
