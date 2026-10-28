"""Extraction of EC numbers and database cross-references from cobrapy reactions."""

from __future__ import annotations

import re
from typing import Any, Iterable

from .ec import normalize_ec

# annotation keys (identifiers.org prefixes) that hold EC numbers
EC_NAMESPACES = ("ec-code", "ec", "ec-number", "brenda")

# annotation key -> canonical namespace of the cross-reference
XREF_NAMESPACES = {
    "bigg.reaction": "bigg.reaction",
    "biggr": "bigg.reaction",  # models from bigg.bio

    "kegg.reaction": "kegg.reaction",
    "seed.reaction": "seed.reaction",
    "seed": "seed.reaction",  # legacy identifiers.org/seed, only rxn IDs are kept
    "rhea": "rhea",
    "metacyc.reaction": "metacyc.reaction",
    "biocyc": "metacyc.reaction",
}

# patterns of valid identifiers per canonical namespace
_ID_PATTERNS = {
    "bigg.reaction": re.compile(r"\S+"),
    "kegg.reaction": re.compile(r"R\d{5}"),
    "seed.reaction": re.compile(r"rxn\d{5}"),
    "rhea": re.compile(r"\d+"),
    "metacyc.reaction": re.compile(r"[^\s:]+"),
}

# bigg.bio gives reactions and metabolites that are not in its universal namespace model specific
# IDs, e.g. __iYL1228__PGI or __iYL1228__g6p_c
_MODEL_PREFIX = re.compile(r"^__.+?__(?=.)")

# COBRA Toolbox style notes, e.g. <p>EC Number: 1.1.1.1</p>
_NOTES_EC_KEY = re.compile(r"ec[\s_-]?numbers?|protein[\s_]class", re.IGNORECASE)
_NOTES_SPLIT = re.compile(r"[\s,;|]+")


def annotation_values(annotation: dict[str, Any], key: str) -> list[str]:
    """Values of an annotation entry as list (cobrapy stores either str or list)."""
    value = annotation.get(key)
    if value is None:
        return []
    if isinstance(value, str):
        return [value]
    return [str(v) for v in value]


def _unique(values: Iterable[str]) -> list[str]:
    return list(dict.fromkeys(values))


def strip_model_prefix(identifier: str) -> str:
    """The original ID of a model specific BiGG ID of bigg.bio: ``__iYL1228__PGI`` gives ``PGI``."""
    return _MODEL_PREFIX.sub("", identifier)


def normalize_xref(namespace: str, raw: str) -> str | None:
    """Strip prefixes like ``RHEA:``, ``rn:``, ``META:`` or ``__iYL1228__`` and validate the identifier."""
    value = raw.strip()
    if namespace == "kegg.reaction":
        value = re.sub(r"^(kegg\.reaction|rn):", "", value, flags=re.IGNORECASE)
    elif namespace == "rhea":
        value = re.sub(r"^rhea:", "", value, flags=re.IGNORECASE)
    elif namespace == "metacyc.reaction":
        # biocyc identifiers carry the PGDB as prefix, e.g. META:PGLUCISOM-RXN
        value = value.rsplit(":", 1)[-1]
    elif value.lower().startswith(namespace + ":"):
        value = value[len(namespace) + 1 :]
    if namespace == "bigg.reaction":
        value = strip_model_prefix(value)
    pattern = _ID_PATTERNS.get(namespace)
    if pattern is None or not pattern.fullmatch(value):
        return None
    return value


def annotated_ecs(reaction) -> list[str]:
    """EC numbers in the annotation of `reaction`."""
    ecs = (
        normalize_ec(v)
        for key in EC_NAMESPACES
        for v in annotation_values(reaction.annotation, key)
    )
    return _unique(ec for ec in ecs if ec)


def is_ec_notes_key(key: str) -> bool:
    """Whether a key of COBRA style notes holds EC numbers (``EC Number``, ``PROTEIN_CLASS``)."""
    return bool(_NOTES_EC_KEY.fullmatch(key.strip()))


def notes_ecs(reaction) -> list[str]:
    """EC numbers in COBRA style notes (``EC Number``, ``PROTEIN_CLASS``) of `reaction`."""
    ecs = (
        normalize_ec(token)
        for key, value in (reaction.notes or {}).items()
        if is_ec_notes_key(key)
        for token in _NOTES_SPLIT.split(str(value))
        if token
    )
    return _unique(ec for ec in ecs if ec)


def xrefs(reaction) -> dict[str, list[str]]:
    """Cross-references in the annotation of `reaction` as {namespace: [identifiers]}."""
    result: dict[str, list[str]] = {}
    for key, namespace in XREF_NAMESPACES.items():
        for raw in annotation_values(reaction.annotation, key):
            identifier = normalize_xref(namespace, raw)
            if identifier is not None:
                result.setdefault(namespace, [])
                if identifier not in result[namespace]:
                    result[namespace].append(identifier)
    return result
