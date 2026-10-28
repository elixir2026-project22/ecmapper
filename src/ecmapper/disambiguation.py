"""Selection of the EC number(s) of a reaction if several were found."""

from __future__ import annotations

from typing import Iterable

from .ec import is_complete

#: values of :attr:`ecmapper.ReactionECMapping.selection`
SELECTIONS = ("annotation", "notes", "single", "bkms", "ambiguous", "none", "boundary")


def covers(general: str, specific: str) -> bool:
    """Whether the (incomplete) EC number `general` includes `specific`, e.g. ``1.1.1.-`` and ``1.1.1.1``."""
    return general != specific and all(g in ("-", s) for g, s in zip(general.split("."), specific.split(".")))


def drop_covered(ecs: Iterable[str]) -> list[str]:
    """Remove incomplete EC numbers that include another of the EC numbers, sorted."""
    ecs = set(ecs)
    return sorted(ec for ec in ecs if is_complete(ec) or not any(covers(ec, other) for other in ecs))


def select(mapping, bkms=None) -> tuple[list[str], str]:
    """Select the EC number(s) of a :class:`ecmapper.ReactionECMapping`.

    Returns the selected EC numbers and how they were selected (see :data:`SELECTIONS`):

    1. ``annotation``/``notes``: the EC numbers that are already in the model are kept.
    2. ``single``: only one EC number was found (incomplete EC numbers that are
       included in a found complete one, e.g. ``1.1.1.-`` and ``1.1.1.1``, do not count).
    3. ``bkms``: all BKMS-react reactions that match the KEGG or MetaCyc
       reactions of the reaction (from the model or found as cross-reference)
       have the same EC number; `bkms` is the :class:`~ecmapper.sources.BKMSSource`.
    4. ``ambiguous``: all EC numbers found are kept.
    """
    if mapping.origin in ("none", "boundary"):
        return [], mapping.origin
    if mapping.origin in ("annotation", "notes"):
        return sorted({e.ec for e in mapping.evidence if e.source == mapping.origin}), mapping.origin
    candidates = drop_covered(mapping.ec_numbers)
    if len(candidates) == 1:
        return candidates, "single"
    if bkms is not None:
        consensus = drop_covered(bkms.consensus(mapping.identifiers))
        if len(consensus) == 1:
            return consensus, "bkms"
    return candidates, "ambiguous"
