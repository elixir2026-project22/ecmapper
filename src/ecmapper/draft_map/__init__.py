"""Draft metabolic maps of a model from the BRENDA map.

Starting from the matches of a model to the BRENDA map (:mod:`ecmapper.brenda_map`), the map is reduced to the
subsystems and enzyme nodes with model reactions, completed with model reactions that connect drawn metabolites,
condensed, and the subsystems are arranged with few crossings of inter-pathway links. Output: SVG in the style of the
BRENDA map, Escher map, CNApy box positions. See DRAFT_MAP.md; command line: :mod:`.prepare`, then :mod:`.run`.
"""
