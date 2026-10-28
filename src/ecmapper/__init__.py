"""Map reactions of constraint-based models to (BRENDA) EC numbers.

Example
-------
>>> import cobra
>>> from ecmapper import ECMapper, to_dataframe
>>> model = cobra.io.read_sbml_model("model.xml")
>>> mappings = ECMapper().map_model(model)
>>> to_dataframe(mappings)
"""

from .disambiguation import SELECTIONS, TIE_BREAKERS, Disambiguator
from .ec import is_complete, normalize_ec
from .evaluation import evaluate, evaluation_summary
from .metabolites import MetaboliteComparison
from .mapper import (
    ECEvidence,
    ECMapper,
    ReactionECMapping,
    annotate_model,
    map_model,
    summary,
    to_dataframe,
)
from .sources import RESOURCES, SOURCES, FileStatus, check_updates, default_sources, update_sources

__version__ = "0.1.0"

__all__ = [
    "RESOURCES",
    "SELECTIONS",
    "SOURCES",
    "TIE_BREAKERS",
    "Disambiguator",
    "MetaboliteComparison",
    "ECEvidence",
    "ECMapper",
    "FileStatus",
    "ReactionECMapping",
    "annotate_model",
    "check_updates",
    "default_sources",
    "evaluate",
    "evaluation_summary",
    "is_complete",
    "map_model",
    "normalize_ec",
    "summary",
    "to_dataframe",
    "update_sources",
]
