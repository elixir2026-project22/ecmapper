"""Resources that map reaction identifiers to EC numbers."""

from __future__ import annotations

from typing import Iterable

from .base import FILE_STATUSES, ECSource, FileStatus, default_cache_dir
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


def _as_sources(sources: Iterable[ECSource | str] | None, **kwargs) -> list[ECSource]:
    if sources is None:
        return default_sources(**kwargs)
    return [s if isinstance(s, ECSource) else SOURCES[s](**kwargs) for s in sources]


def check_updates(sources: Iterable[ECSource | str] | None = None, **kwargs) -> list[FileStatus]:
    """Check whether newer versions of the cached files of `sources` (default: all) are available.

    `sources` are instances or names, `kwargs` (e.g. `cache_dir`) are passed to sources
    created from names. Nothing in the cache is changed.
    """
    return [status for source in _as_sources(sources, **kwargs) for status in source.check_updates()]


def update_sources(
    sources: Iterable[ECSource | str] | None = None, force: bool = False, **kwargs
) -> list[FileStatus]:
    """Download newer versions of the cached files of `sources` (default: all).

    Only files that changed (or are missing) are downloaded, with `force` all files.
    """
    return [status for source in _as_sources(sources, **kwargs) for status in source.update(force=force)]


__all__ = [
    "FILE_STATUSES",
    "SOURCES",
    "BKMSSource",
    "BiGGrSource",
    "ECSource",
    "FileStatus",
    "KEGGSource",
    "ModelSEEDSource",
    "RheaSource",
    "check_updates",
    "default_cache_dir",
    "default_sources",
    "update_sources",
]
