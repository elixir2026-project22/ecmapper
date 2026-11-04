"""GUI for the assignment of the subsystems of a draft map to maps (--split) and their arrangement.

usage: python -m ecmapper.draft_map.gui PREP.pkl BRENDA_map_reactions.json BRENDA_map.svg OUT_DIR [--assignment XLSX]
[options of ecmapper.draft_map.run]

:mod:`.session` holds the state (no Qt), :mod:`.items` draws the maps, :mod:`.window` is the main window.
"""
