"""Mapping of model reactions to EC numbers."""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from typing import Iterable, Sequence

import pandas as pd

from . import annotation
from .disambiguation import SELECTIONS, Disambiguator
from .ec import is_complete
from .metabolites import MetaboliteComparison
from .sources import (
    BiGGrMetabolites,
    BKMSSource,
    CachedResource,
    ECSource,
    EnzymeNomenclature,
    FileStatus,
    KEGGReactionCompounds,
    KEGGSource,
    default_sources,
)

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
    identifiers : dict
        All identifiers of the reaction that were used, {namespace: [identifiers]}:
        the cross-references in the model and those found as cross-reference
        in a resource (only filled if the resources were queried).
    selected : list of str
        The EC numbers selected from all found ones, see :class:`ecmapper.disambiguation.Disambiguator`.
    selection : str
        How they were selected (see :data:`ecmapper.disambiguation.SELECTIONS`).
    scores : dict
        Scores of the comparison of metabolites {EC number: score}, None for EC
        numbers without KEGG reactions; empty if no comparison was made.
    """

    reaction_id: str
    origin: str
    evidence: list[ECEvidence] = field(default_factory=list)
    xrefs: dict[str, list[str]] = field(default_factory=dict)
    identifiers: dict[str, list[str]] = field(default_factory=dict)
    selected: list[str] = field(default_factory=list)
    selection: str = ""
    scores: dict[str, float | None] = field(default_factory=dict)

    @property
    def ec_numbers(self) -> list[str]:
        """All EC numbers, sorted."""
        return sorted({e.ec for e in self.evidence})

    @property
    def ambiguous(self) -> bool:
        """True if several EC numbers were found and none could be selected."""
        return self.selection == "ambiguous"

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
    If several EC numbers were found, the ones to use are selected, see
    :class:`ecmapper.disambiguation.Disambiguator`.

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
    bkms_consensus : bool
        Select the EC number on which all matching BKMS-react reactions agree
        (requires the ``bkms`` source).
    current_ec : bool
        Replace transferred EC numbers by their successors and drop deleted ones
        when selecting (ExPASy ENZYME).
    compare_metabolites : bool
        Select EC numbers by comparing the metabolites of a reaction with the
        KEGG reactions of the EC numbers (KEGG and BiGGr metabolites).
    min_score, min_margin : float
        Thresholds of the comparison of metabolites.
    tie_breakers : sequence of str
        ``sources`` and/or ``organism``, applied in this order to the EC numbers
        that remain; none by default.
    organism : str, optional
        Swiss-Prot organism mnemonic (e.g. ``ECOLI``) for the ``organism`` tie-breaker.
    cache_dir, refresh :
        Passed to the default sources and the reference data, see
        :class:`ecmapper.sources.CachedResource`. If sources are given as
        instances and `cache_dir` is not set, the reference data is cached
        where the first source is cached.
    """

    def __init__(
        self,
        sources: Iterable[ECSource | str] | None = None,
        *,
        use_notes: bool = True,
        always_query: bool = False,
        skip_boundary: bool = True,
        follow_xrefs: bool = True,
        bkms_consensus: bool = True,
        current_ec: bool = True,
        compare_metabolites: bool = True,
        min_score: float = 0.6,
        min_margin: float = 0.1,
        tie_breakers: Sequence[str] = (),
        organism: str | None = None,
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
            if cache_dir is None and sources:
                cache_dir = sources[0].cache_dir
        self.sources: dict[str, ECSource] = {}
        for source in sources:
            if source.namespace in self.sources:
                raise ValueError(f"more than one source for namespace {source.namespace}")
            self.sources[source.namespace] = source
        self.use_notes = use_notes
        self.always_query = always_query
        self.skip_boundary = skip_boundary
        self.follow_xrefs = follow_xrefs

        def reference(cls):
            # a source of the mapper is shared with the disambiguation, others are created
            existing = self.sources.get(getattr(cls, "namespace", None))
            return existing if isinstance(existing, cls) else cls(cache_dir=cache_dir, refresh=refresh)

        enzyme = reference(EnzymeNomenclature) if current_ec or "organism" in tie_breakers else None
        comparison = None
        if compare_metabolites:
            comparison = MetaboliteComparison(
                reference(KEGGSource),
                reference(KEGGReactionCompounds),
                reference(BiGGrMetabolites),
                enzyme,
            )
        bkms = next((s for s in self.sources.values() if isinstance(s, BKMSSource)), None)
        self.disambiguator = Disambiguator(
            enzyme=enzyme,
            bkms=bkms if bkms_consensus else None,
            comparison=comparison,
            min_score=min_score,
            min_margin=min_margin,
            tie_breakers=tie_breakers,
            organism=organism,
        )

    @property
    def resources(self) -> list[CachedResource]:
        """All resources with cached files: the sources and the reference data of the disambiguation."""
        resources = list(self.sources.values())
        d = self.disambiguator
        for resource in [d.enzyme, *(d.comparison.resources if d.comparison else [])]:
            if resource is not None and all(resource is not r for r in resources):
                resources.append(resource)
        return resources

    def check_updates(self) -> list[FileStatus]:
        """Check whether newer versions of the resource files are available, see :meth:`CachedResource.check_updates`."""
        return [status for resource in self.resources for status in resource.check_updates()]

    def update(self, force: bool = False) -> list[FileStatus]:
        """Download newer versions of the resource files, see :meth:`CachedResource.update`."""
        statuses = [status for resource in self.resources for status in resource.update(force=force)]
        if self.disambiguator.comparison is not None:
            self.disambiguator.comparison.reset()
        return statuses

    def map_reaction(self, reaction, boundary: bool = False) -> ReactionECMapping:
        """Map one cobrapy reaction, `boundary` marks it as boundary reaction."""
        xrefs = annotation.xrefs(reaction)
        if boundary and self.skip_boundary:
            return ReactionECMapping(reaction.id, "boundary", xrefs=xrefs, selection="boundary")
        origin = "none"
        identifiers = {}
        evidence = [ECEvidence(ec, "annotation") for ec in annotation.annotated_ecs(reaction)]
        if evidence:
            origin = "annotation"
        elif self.use_notes:
            evidence = [ECEvidence(ec, "notes") for ec in annotation.notes_ecs(reaction)]
            if evidence:
                origin = "notes"
        if not evidence or self.always_query:
            identifiers = self.expand_xrefs(xrefs)
            found = self._lookup(identifiers)
            evidence += found
            if found and origin == "none":
                origin = "resources"
        mapping = ReactionECMapping(
            reaction.id, origin, evidence, xrefs, {ns: list(ids) for ns, ids in identifiers.items()}
        )
        selection = self.disambiguator.select(mapping, reaction)
        mapping.selected, mapping.selection, mapping.scores = selection.selected, selection.selection, selection.scores
        return mapping

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

    def expand_xrefs(self, xrefs: dict[str, list[str]]) -> dict[str, dict[str, str]]:
        """Add the cross-references the resources list for `xrefs` {namespace: [identifiers]}.

        Returns {namespace: {identifier: via}}, where `via` is empty for the
        identifiers in `xrefs` and ``namespace:identifier`` of the identifier
        whose cross-reference it is otherwise. Cross-references are only
        followed if :attr:`follow_xrefs` is set.
        """
        result = {ns: dict.fromkeys(ids, "") for ns, ids in xrefs.items()}
        if not self.follow_xrefs:
            return result
        for namespace, identifiers in xrefs.items():
            source = self.sources.get(namespace)
            if source is None:
                continue
            for identifier in identifiers:
                for xref_namespace, xref_ids in sorted(source.xrefs(identifier).items()):
                    for xref_id in sorted(xref_ids):
                        result.setdefault(xref_namespace, {}).setdefault(xref_id, f"{namespace}:{identifier}")
        return result

    def query(self, xrefs: dict[str, list[str]]) -> list[ECEvidence]:
        """Look up the cross-references {namespace: [identifiers]} in the resources."""
        return self._lookup(self.expand_xrefs(xrefs))

    def _lookup(self, identifiers: dict[str, dict[str, str]]) -> list[ECEvidence]:
        evidence = []
        for namespace, ids in identifiers.items():
            source = self.sources.get(namespace)
            if source is None:
                continue
            for identifier, via in ids.items():
                evidence += [ECEvidence(ec, source.name, identifier, via) for ec in sorted(source.lookup(identifier))]
        return evidence


def map_model(model, **kwargs) -> dict[str, ReactionECMapping]:
    """Shortcut for ``ECMapper(**kwargs).map_model(model)``."""
    return ECMapper(**kwargs).map_model(model)


def _join(values: Iterable[str]) -> str:
    return ";".join(values)


def to_dataframe(mappings: dict[str, ReactionECMapping] | Iterable[ReactionECMapping], long: bool = False) -> pd.DataFrame:
    """Tabulate mappings.

    The default (wide) table has one row per reaction with the selected and
    all found EC numbers separated by ``;``. The long table has one row per
    reaction, EC number and evidence, reactions without EC number get one row
    with empty EC.
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
                        "selection": m.selection,
                        "ec": e.ec,
                        "selected": e.ec in m.selected,
                        "score": m.scores.get(e.ec),
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
                    "selection": m.selection,
                    "selected": _join(m.selected),
                    "ec_numbers": _join(m.ec_numbers),
                    "n_ec": len(m.ec_numbers),
                    "sources": _join(m.sources),
                    "scores": _join(f"{ec}:{s:.2f}" for ec, s in m.scores.items() if s is not None),
                    "xrefs": _join(f"{ns}:{i}" for ns, ids in m.xrefs.items() for i in ids),
                }
            )
    columns = (
        ["reaction_id", "origin", "selection", "ec", "selected", "score", "complete", "source", "identifier", "via"]
        if long
        else ["reaction_id", "origin", "selection", "selected", "ec_numbers", "n_ec", "sources", "scores", "xrefs"]
    )
    return pd.DataFrame(rows, columns=columns)


def summary(mappings: dict[str, ReactionECMapping]) -> dict:
    """Counts of reactions per origin, per selection, of ambiguous reactions and per source."""
    values = list(mappings.values())
    origins = Counter(m.origin for m in values)
    selections = Counter(m.selection for m in values)
    per_source = Counter(s for m in values for s in m.sources)
    return {
        "reactions": len(values),
        "origin": {o: origins.get(o, 0) for o in ORIGINS},
        "selection": {s: selections.get(s, 0) for s in SELECTIONS},
        "with_ec": sum(1 for m in values if m.evidence),
        "ambiguous": sum(1 for m in values if m.ambiguous),
        "reactions_per_source": dict(per_source.most_common()),
    }


def annotate_model(model, mappings: dict[str, ReactionECMapping], overwrite: bool = False) -> int:
    """Write the selected EC numbers found in the resources to ``reaction.annotation["ec-code"]``.

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
        reaction.annotation["ec-code"] = m.selected
        changed += 1
    return changed
