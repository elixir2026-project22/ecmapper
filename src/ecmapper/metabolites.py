"""Comparison of the metabolites of a model reaction with the reactions of EC numbers."""

from __future__ import annotations

import re
from typing import Iterable

from .annotation import annotation_values, strip_model_prefix
from .sources import BiGGrMetabolites, EnzymeNomenclature, KEGGReactionCompounds, KEGGSource

PROTON = "C00080"

#: KEGG compounds that are cofactors or otherwise take part in many reactions
CURRENCY = frozenset(
    {
        "C00001",  # H2O
        "C00002",  # ATP
        "C00008",  # ADP
        "C00020",  # AMP
        "C00003",  # NAD+
        "C00004",  # NADH
        "C00005",  # NADPH
        "C00006",  # NADP+
        "C00009",  # orthophosphate
        "C00013",  # diphosphate
        "C00010",  # CoA
        "C00011",  # CO2
        "C00014",  # NH3
        "C01342",  # NH4+
        "C00007",  # O2
        "C00027",  # H2O2
        "C00044",  # GTP
        "C00035",  # GDP
        "C00144",  # GMP
        "C00063",  # CTP
        "C00112",  # CDP
        "C00055",  # CMP
        "C00075",  # UTP
        "C00015",  # UDP
        "C00105",  # UMP
        "C00016",  # FAD
        "C01352",  # FADH2
        "C00399",  # ubiquinone
        "C00390",  # ubiquinol
    }
)

_KEGG_COMPOUND = re.compile(r"C\d{5}")
_COMPARTMENT_SUFFIX = re.compile(r"_[a-z][a-z0-9]?$")


class MetaboliteComparison:
    """Scores how well the metabolites of a model reaction match the KEGG reactions of EC numbers.

    The metabolites of the model reaction are mapped to KEGG compounds by their
    ``kegg.compound`` annotation and, via their ``bigg.metabolite`` annotation
    or their ID without compartment, by the BiGGr metabolites download. Protons
    and metabolites without KEGG compound are ignored.

    The score of an EC number is the highest weighted Jaccard index between
    the metabolites and the compounds of one of the KEGG reactions of the EC
    number: a metabolite matches if one of its KEGG compounds is in the
    reaction, currency metabolites (:data:`CURRENCY`) count `currency_weight`,
    all others 1. The sides of the reactions are not compared.

    Parameters
    ----------
    kegg : KEGGSource
        For the KEGG reactions of the EC numbers.
    compounds : KEGGReactionCompounds
    bigg_metabolites : BiGGrMetabolites, optional
    enzyme : EnzymeNomenclature, optional
        If given, the KEGG reactions of transferred EC numbers count for the
        EC numbers they were transferred to.
    currency_weight : float
    """

    def __init__(
        self,
        kegg: KEGGSource,
        compounds: KEGGReactionCompounds,
        bigg_metabolites: BiGGrMetabolites | None = None,
        enzyme: EnzymeNomenclature | None = None,
        currency_weight: float = 0.5,
    ):
        self.kegg = kegg
        self.compounds = compounds
        self.bigg_metabolites = bigg_metabolites
        self.enzyme = enzyme
        self.currency_weight = currency_weight
        self._ec_reactions: dict[str, set[str]] | None = None

    @property
    def resources(self) -> list:
        return [r for r in (self.kegg, self.compounds, self.bigg_metabolites, self.enzyme) if r is not None]

    def reset(self) -> None:
        """Forget the index of KEGG reactions per EC number (after an update of the resources)."""
        self._ec_reactions = None

    def reference_reactions(self, ec: str) -> set[str]:
        """The KEGG reactions of the EC number `ec`."""
        if self._ec_reactions is None:
            index: dict[str, set[str]] = {}
            for reaction, ecs in self.kegg.table.items():
                for e in ecs:
                    for current in self.enzyme.current(e) if self.enzyme else {e}:
                        index.setdefault(current, set()).add(reaction)
            self._ec_reactions = index
        return set(self._ec_reactions.get(ec, ()))

    def metabolite_compounds(self, metabolite) -> set[str]:
        """KEGG compound IDs of a cobrapy metabolite."""
        compounds = {
            v.removeprefix("kegg.compound:").removeprefix("cpd:")
            for v in annotation_values(metabolite.annotation, "kegg.compound")
        }
        if self.bigg_metabolites is not None:
            for bigg_id in self.bigg_ids(metabolite):
                compounds |= self.bigg_metabolites.xrefs(bigg_id).get("kegg.compound", set())
        return {c for c in compounds if _KEGG_COMPOUND.fullmatch(c) and c != PROTON}

    @staticmethod
    def bigg_ids(metabolite) -> set[str]:
        """Possible universal BiGG IDs of a cobrapy metabolite: annotation and ID without compartment.

        Models from bigg.bio annotate compartmentalised metabolites with charge
        (``biggr: g6p_c:-2``), which gives ``g6p``, and give metabolites that
        are not in the universal namespace model specific IDs
        (``__iYL1228__g6p_c``), which give the original ID.
        """
        bigg_ids = set(annotation_values(metabolite.annotation, "bigg.metabolite"))
        suffix = f"_{metabolite.compartment}" if metabolite.compartment else None
        identifiers = [metabolite.id, *(v.split(":")[0] for v in annotation_values(metabolite.annotation, "biggr"))]
        for identifier in map(strip_model_prefix, identifiers):
            if suffix and identifier.endswith(suffix):
                bigg_ids.add(identifier[: -len(suffix)])
            else:
                # compartments with other IDs than the suffixes of the metabolites (Cytosol, glc__D_c)
                bigg_ids.add(_COMPARTMENT_SUFFIX.sub("", identifier))
        return bigg_ids

    def reaction_compounds(self, reaction) -> list[frozenset[str]]:
        """KEGG compounds of the metabolites of a cobrapy reaction, one set per metabolite (without duplicates)."""
        result = []
        for metabolite in reaction.metabolites:
            compounds = frozenset(self.metabolite_compounds(metabolite))
            if compounds and compounds not in result:
                result.append(compounds)
        return result

    def score(self, metabolites: list[frozenset[str]], reference: Iterable[str]) -> float:
        """Weighted Jaccard index between metabolites (sets of KEGG compounds) and the compounds of a reaction."""
        reference = set(reference) - {PROTON}
        matched, weight_matched, weight_union = set(), 0.0, 0.0
        for compounds in metabolites:
            weight = self.currency_weight if compounds & CURRENCY else 1.0
            weight_union += weight
            hits = compounds & reference
            if hits:
                weight_matched += weight
                matched |= hits
        for compound in reference - matched:
            weight_union += self.currency_weight if compound in CURRENCY else 1.0
        return weight_matched / weight_union if weight_union else 0.0

    def scores(self, reaction, ecs: Iterable[str]) -> dict[str, float | None]:
        """Score of each EC number for the cobrapy `reaction`.

        None for EC numbers without KEGG reactions. Empty if no metabolite of the
        reaction could be mapped to a KEGG compound.
        """
        metabolites = self.reaction_compounds(reaction)
        if not metabolites:
            return {}
        result = {}
        for ec in ecs:
            references = self.reference_reactions(ec)
            result[ec] = (
                max(self.score(metabolites, self.compounds.compounds(r)) for r in references) if references else None
            )
        return result
