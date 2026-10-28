"""Resources that map reaction identifiers to EC numbers, and reference data for the disambiguation."""

from __future__ import annotations

from typing import Iterable

from .base import FILE_STATUSES, CachedResource, ECSource, FileStatus, default_cache_dir
from .biggr import BiGGrMetabolites, BiGGrSource
from .bkms import BKMSSource
from .enzyme import EnzymeEntry, EnzymeNomenclature
from .kegg import KEGGReactionCompounds, KEGGSource
from .rhea import RheaSource
from .seed import ModelSEEDSource

#: all available EC sources by name
SOURCES: dict[str, type[ECSource]] = {
    cls.name: cls for cls in (KEGGSource, ModelSEEDSource, RheaSource, BKMSSource, BiGGrSource)
}

#: reference data used to select EC numbers, by name
REFERENCES: dict[str, type[CachedResource]] = {
    cls.name: cls for cls in (EnzymeNomenclature, KEGGReactionCompounds, BiGGrMetabolites)
}

#: all resources with cached files, by name
RESOURCES: dict[str, type[CachedResource]] = {**SOURCES, **REFERENCES}


def default_sources(names=None, **kwargs) -> list[ECSource]:
    """Instantiate the sources `names` (default: all), `kwargs` are passed to each source."""
    names = list(SOURCES) if names is None else list(names)
    unknown = set(names) - set(SOURCES)
    if unknown:
        raise ValueError(f"unknown source(s) {sorted(unknown)}, available are {list(SOURCES)}")
    return [SOURCES[name](**kwargs) for name in names]


def _as_resources(resources: Iterable[CachedResource | str] | None, **kwargs) -> list[CachedResource]:
    if resources is None:
        resources = list(RESOURCES)
    resources = list(resources)
    unknown = {r for r in resources if isinstance(r, str)} - set(RESOURCES)
    if unknown:
        raise ValueError(f"unknown resource(s) {sorted(unknown)}, available are {list(RESOURCES)}")
    return [r if isinstance(r, CachedResource) else RESOURCES[r](**kwargs) for r in resources]


def check_updates(resources: Iterable[CachedResource | str] | None = None, **kwargs) -> list[FileStatus]:
    """Check whether newer versions of the cached files of `resources` (default: all) are available.

    `resources` are instances or names (of :data:`RESOURCES`), `kwargs` (e.g.
    `cache_dir`) are passed to resources created from names. Nothing in the
    cache is changed.
    """
    return [status for resource in _as_resources(resources, **kwargs) for status in resource.check_updates()]


def update_sources(
    resources: Iterable[CachedResource | str] | None = None, force: bool = False, **kwargs
) -> list[FileStatus]:
    """Download newer versions of the cached files of `resources` (default: all).

    Only files that changed (or are missing) are downloaded, with `force` all files.
    """
    return [status for resource in _as_resources(resources, **kwargs) for status in resource.update(force=force)]


__all__ = [
    "FILE_STATUSES",
    "REFERENCES",
    "RESOURCES",
    "SOURCES",
    "BKMSSource",
    "BiGGrMetabolites",
    "BiGGrSource",
    "CachedResource",
    "ECSource",
    "EnzymeEntry",
    "EnzymeNomenclature",
    "FileStatus",
    "KEGGReactionCompounds",
    "KEGGSource",
    "ModelSEEDSource",
    "RheaSource",
    "check_updates",
    "default_cache_dir",
    "default_sources",
    "update_sources",
]
