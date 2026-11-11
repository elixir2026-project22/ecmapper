"""PWL (Pathway Layout Language) 2.0: text descriptions of metabolic subsystems, see ``PWL.md``.

- :func:`parse` / :func:`parse_file` read a description (no dependencies),
- :func:`detect_motifs` finds motifs in it (cycles, substrate cycles, shunts, ...),
- :func:`ecmapper.pwl.brenda.generate` writes descriptions of the subsystems of a BRENDA map
  export (needs numpy and networkx), :func:`ecmapper.pwl.check.check` validates them.
"""

from .motifs import detect_motifs, format_motif
from .syntax import (
    VERSION,
    Compound,
    Document,
    Motif,
    Participant,
    PWLError,
    Reaction,
    parse,
    parse_file,
    validate,
)

__all__ = [
    "VERSION",
    "Compound",
    "Document",
    "Motif",
    "PWLError",
    "Participant",
    "Reaction",
    "detect_motifs",
    "format_motif",
    "parse",
    "parse_file",
    "validate",
]
