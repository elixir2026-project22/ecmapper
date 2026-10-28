"""Resources that map reaction identifiers to EC numbers."""

from __future__ import annotations

from .base import ECSource, default_cache_dir
from .biggr import BiGGrSource
from .bkms import BKMSSource
from .kegg import KEGGSource
from .rhea import RheaSource
from .seed import ModelSEEDSource

#: all available sources by name
SOURCES: dict[str, type[ECSource]] = {
    cls.name: cls for cls in (KEGGSource, ModelSEEDSource, RheaSource, BKMSSource, BiGGrSource)
}


def default_sources(names=None, **kwargs) -> list[ECSource]:
    """Instantiate the sources `names` (default: all), `kwargs` are passed to each source."""
    names = list(SOURCES) if names is None else list(names)
    unknown = set(names) - set(SOURCES)
    if unknown:
        raise ValueError(f"unknown source(s) {sorted(unknown)}, available are {list(SOURCES)}")
    return [SOURCES[name](**kwargs) for name in names]


__all__ = [
    "SOURCES",
    "BKMSSource",
    "BiGGrSource",
    "ECSource",
    "KEGGSource",
    "ModelSEEDSource",
    "RheaSource",
    "default_cache_dir",
    "default_sources",
]
