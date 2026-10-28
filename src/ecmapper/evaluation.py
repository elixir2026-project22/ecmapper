"""Evaluation of the mapping with the EC numbers that a model already has."""

from __future__ import annotations

import pandas as pd

from . import annotation
from .disambiguation import SELECTIONS, covers
from .mapper import ECMapper

#: values of the column ``outcome`` of :func:`evaluate`
OUTCOMES = ("consistent", "overlap", "disjoint", "not found")


def _strip(model, bigg_only: bool) -> None:
    """Remove the EC numbers from `model` (in place); with `bigg_only` keep only BiGG annotations."""
    for reaction in model.reactions:
        if bigg_only:
            bigg = reaction.annotation.get("bigg.reaction", reaction.id)
            reaction.annotation = {"bigg.reaction": bigg}
        else:
            for key in annotation.EC_NAMESPACES:
                reaction.annotation.pop(key, None)
        reaction.notes = {k: v for k, v in (reaction.notes or {}).items() if not annotation.is_ec_notes_key(k)}
    if bigg_only:
        for metabolite in model.metabolites:
            metabolite.annotation = {k: v for k, v in metabolite.annotation.items() if k == "bigg.metabolite"}


def evaluate(model, mapper: ECMapper | None = None, bigg_only: bool = False) -> pd.DataFrame:
    """Map a model whose EC numbers were removed and compare the result with them.

    The EC numbers in the annotation and notes of the reactions are the
    reference (made current if the mapper replaces transferred EC numbers). A
    copy of the model without them is mapped by `mapper` (default:
    ``ECMapper()``); with `bigg_only`, the copy keeps only the BiGG annotations
    of reactions and metabolites. Returns one row per reaction with reference
    EC numbers, ``outcome`` is ``consistent`` if all selected EC numbers are in
    the reference, ``overlap`` if some are, ``disjoint`` if none is and
    ``not found`` if no EC number was found. A selected EC number that is
    included in an incomplete reference EC number (``2.5.1.141`` in
    ``2.5.1.-``) counts as in the reference.
    """
    mapper = mapper or ECMapper()
    reference = {}
    for reaction in model.reactions:
        ecs = annotation.annotated_ecs(reaction) + annotation.notes_ecs(reaction)
        if ecs:
            reference[reaction.id] = set(mapper.disambiguator.current(ecs) or ecs)
    stripped = model.copy()
    _strip(stripped, bigg_only)
    mappings = mapper.map_model(stripped, [r for r in stripped.reactions if r.id in reference])

    rows = []
    for reaction_id, ecs in reference.items():
        m = mappings[reaction_id]
        matching = [ec for ec in m.selected if ec in ecs or any(covers(r, ec) for r in ecs)]
        if not m.selected:
            outcome = "not found"
        elif len(matching) == len(m.selected):
            outcome = "consistent"
        elif matching:
            outcome = "overlap"
        else:
            outcome = "disjoint"
        rows.append(
            {
                "reaction_id": reaction_id,
                "reference": ";".join(sorted(ecs)),
                "selection": m.selection,
                "selected": ";".join(m.selected),
                "candidates": ";".join(m.ec_numbers),
                "outcome": outcome,
            }
        )
    return pd.DataFrame(rows, columns=["reaction_id", "reference", "selection", "selected", "candidates", "outcome"])


def evaluation_summary(table: pd.DataFrame) -> pd.DataFrame:
    """Number of reactions per selection (rows) and outcome (columns)."""
    summary = pd.crosstab(table["selection"], table["outcome"])
    summary = summary.reindex(
        index=[s for s in SELECTIONS if s in summary.index], columns=list(OUTCOMES), fill_value=0
    )
    summary["total"] = summary.sum(axis=1)
    summary.loc["total"] = summary.sum()
    return summary
