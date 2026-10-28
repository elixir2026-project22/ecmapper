"""Validation and normalisation of EC numbers."""

from __future__ import annotations

import re

_PREFIX = re.compile(r"^EC[\s:_-]*", re.IGNORECASE)
_SERIAL = re.compile(r"[nB]?\d+")


def normalize_ec(raw: str) -> str | None:
    """Return `raw` as canonical four level EC number or None if it is not one.

    Accepted are official EC numbers, preliminary IUBMB numbers (``1.1.1.n1``),
    BRENDA's provisional numbers (``2.3.1.B43``) and incomplete EC numbers.
    Incomplete numbers are padded with ``-`` (``2.4.1`` becomes ``2.4.1.-``).
    Prefixes like ``EC:`` and surrounding whitespace are removed. Everything
    else, e.g. transporter classification numbers (``3.A.1.7.-``), ``SPONTANEOUS``
    or the letter suffixes used by some databases (``7.3.2.a``) is rejected.
    """
    s = _PREFIX.sub("", raw.strip()).strip()
    parts = s.split(".")
    if not 2 <= len(parts) <= 4:
        return None
    parts += ["-"] * (4 - len(parts))
    if not parts[0].isdigit():
        return None
    incomplete = False
    for level, part in enumerate(parts[1:], start=2):
        if part == "-":
            incomplete = True
        elif incomplete:
            return None
        elif level < 4 and not part.isdigit():
            return None
        elif level == 4 and not _SERIAL.fullmatch(part):
            return None
    return ".".join(parts)


def is_complete(ec: str) -> bool:
    """Whether the (normalised) EC number specifies all four levels."""
    return "-" not in ec
