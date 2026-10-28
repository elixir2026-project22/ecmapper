"""Map reactions of constraint-based models to (BRENDA) EC numbers.

Example
-------
>>> import cobra
>>> from ecmapper import ECMapper, to_dataframe
>>> model = cobra.io.read_sbml_model("model.xml")
>>> mappings = ECMapper().map_model(model)
>>> to_dataframe(mappings)
"""

from .ec import is_complete, normalize_ec
from .mapper import (
    ECEvidence,
    ECMapper,
    ReactionECMapping,
    annotate_model,
    map_model,
    summary,
    to_dataframe,
)
from .sources import SOURCES, default_sources

__version__ = "0.1.0"

__all__ = [
    "SOURCES",
    "ECEvidence",
    "ECMapper",
    "ReactionECMapping",
    "annotate_model",
    "default_sources",
    "is_complete",
    "map_model",
    "normalize_ec",
    "summary",
    "to_dataframe",
]
