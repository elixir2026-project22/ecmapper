"""Selection of the EC number(s) of a reaction if several were found."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Sequence

from .ec import is_complete

#: values of :attr:`ecmapper.ReactionECMapping.selection`
SELECTIONS = (
    "annotation",
    "notes",
    "single",
    "bkms",
    "metabolites",
    "sources",
    "organism",
    "ambiguous",
    "none",
    "boundary",
)

#: tie-breakers that can be enabled in :class:`Disambiguator`
TIE_BREAKERS = ("sources", "organism")


def covers(general: str, specific: str) -> bool:
    """Whether the (incomplete) EC number `general` includes `specific`, e.g. ``1.1.1.-`` and ``1.1.1.1``."""
    return general != specific and all(g in ("-", s) for g, s in zip(general.split("."), specific.split(".")))


def drop_covered(ecs: Iterable[str]) -> list[str]:
    """Remove incomplete EC numbers that include another of the EC numbers, sorted."""
    ecs = set(ecs)
    return sorted(ec for ec in ecs if is_complete(ec) or not any(covers(ec, other) for other in ecs))


@dataclass
class Selection:
    """Result of :meth:`Disambiguator.select`."""

    selected: list[str]
    selection: str
    scores: dict[str, float | None] = field(default_factory=dict)


@dataclass
class Disambiguator:
    """Selects the EC number(s) of a :class:`ecmapper.ReactionECMapping`.

    The rules are applied in this order:

    0. Transferred EC numbers are replaced by their successors, deleted ones
       dropped (if `enzyme` is given), incomplete EC numbers that include a
       found complete one (``2.7.1.-`` and ``2.7.1.165``) are dropped.
    1. ``annotation``/``notes``: the EC numbers that are already in the model are kept.
    2. ``single``: only one EC number was found.
    3. ``bkms``: all BKMS-react reactions that match the KEGG or MetaCyc
       reactions of the reaction have the same EC number (if `bkms` is given).
    4. ``metabolites``: the metabolites of the reaction match the KEGG
       reactions of one EC number clearly best (if `comparison` is given): its
       score is at least `min_score` and more than `min_margin` higher than
       that of all others. If several EC numbers are within `min_margin` of the
       best score, only these are kept.
    5. Tie-breakers, in the order of `tie_breakers`, each keeps the best of
       the remaining EC numbers and selects it if it is the only one:
       ``sources``: the EC numbers supported by most sources;
       ``organism``: the EC numbers of which ENZYME lists a Swiss-Prot protein
       of `organism` (a Swiss-Prot organism mnemonic like ``ECOLI``).
    6. ``ambiguous``: the remaining EC numbers are kept.

    Parameters
    ----------
    enzyme : EnzymeNomenclature, optional
    bkms : BKMSSource, optional
    comparison : MetaboliteComparison, optional
    min_score, min_margin : float
    tie_breakers : sequence of str
        Values of :data:`TIE_BREAKERS`, none by default.
    organism : str, optional
        Required for the ``organism`` tie-breaker.
    """

    enzyme: object = None
    bkms: object = None
    comparison: object = None
    min_score: float = 0.6
    min_margin: float = 0.1
    tie_breakers: Sequence[str] = ()
    organism: str | None = None

    def __post_init__(self):
        self.tie_breakers = tuple(self.tie_breakers)
        unknown = set(self.tie_breakers) - set(TIE_BREAKERS)
        if unknown:
            raise ValueError(f"unknown tie-breaker(s) {sorted(unknown)}, available are {list(TIE_BREAKERS)}")
        if "organism" in self.tie_breakers and not (self.organism and self.enzyme is not None):
            raise ValueError("the organism tie-breaker requires an organism and the ENZYME data")

    def current(self, ecs: Iterable[str]) -> list[str]:
        """Step 0: current EC numbers without the incomplete ones that include another one."""
        ecs = set(ecs)
        if self.enzyme is not None:
            ecs = set().union(*(self.enzyme.current(ec) for ec in ecs)) if ecs else set()
        return drop_covered(ecs)

    def select(self, mapping, reaction=None) -> Selection:
        """Select the EC numbers of `mapping`; `reaction` (cobrapy) is needed for the comparison of metabolites."""
        if mapping.origin in ("none", "boundary"):
            return Selection([], mapping.origin)
        if mapping.origin in ("annotation", "notes"):
            ecs = {e.ec for e in mapping.evidence if e.source == mapping.origin}
            return Selection(self.current(ecs) or sorted(ecs), mapping.origin)

        candidates = self.current(mapping.ec_numbers) or drop_covered(mapping.ec_numbers)
        if len(candidates) == 1:
            return Selection(candidates, "single")
        if self.bkms is not None:
            consensus = self.current(self.bkms.consensus(mapping.identifiers))
            if len(consensus) == 1:
                return Selection(consensus, "bkms")

        scores = {}
        if self.comparison is not None and reaction is not None:
            scores = self.comparison.scores(reaction, candidates)
            known = {ec: s for ec, s in scores.items() if s is not None}
            best = max(known.values(), default=0.0)
            if best >= self.min_score:
                candidates = sorted(ec for ec, s in known.items() if s >= best - self.min_margin)
                if len(candidates) == 1:
                    return Selection(candidates, "metabolites", scores)

        for tie_breaker in self.tie_breakers:
            if tie_breaker == "sources":
                support = {ec: self._support(mapping, ec) for ec in candidates}
                most = max(support.values())
                candidates = [ec for ec in candidates if support[ec] == most]
            elif tie_breaker == "organism":
                candidates = [ec for ec in candidates if self.organism in self.enzyme.organisms(ec)] or candidates
            if len(candidates) == 1:
                return Selection(candidates, tie_breaker, scores)
        return Selection(candidates, "ambiguous", scores)

    def _support(self, mapping, ec: str) -> int:
        """Number of sources whose EC numbers (made current) include `ec`."""
        return len({e.source for e in mapping.evidence if ec in (self.current([e.ec]) or [e.ec])})


def select(mapping, bkms=None) -> tuple[list[str], str]:
    """Select the EC numbers of `mapping` with the rules 1 to 3 of :class:`Disambiguator`."""
    result = Disambiguator(bkms=bkms).select(mapping)
    return result.selected, result.selection
