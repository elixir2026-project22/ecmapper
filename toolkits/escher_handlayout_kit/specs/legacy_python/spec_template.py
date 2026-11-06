"""TEMPLATE for a new hand-layout spec.  Copy to specs/spec_<mapname>.py and fill in.

The pipeline calls build(L) with an escher_handlayout.engine.Layout already loaded from the input map.
You must (1) register the subsystems, (2) optionally split/merge nodes, (3) give coordinates to every PRIMARY
metabolite node (L.setp), (4) call L.reaction(...) for EVERY reaction, (5) set a title position per subsystem.
Anything missing is reported by `python run.py build ...` in one list (SPEC ERROR) so you can fix it in one pass.

Coordinates: x to the right, y downwards, Escher map units. Primary node radius 20, secondary 10.
Typical spacing: 300 between consecutive nodes of a chain; cofactor dots 110-150 away from the reaction line;
lanes of parallel reactions 110-150 apart.  See docs/GUIDE.md for recipes and the full API.
Ids are the string node/reaction ids of the INPUT map (see analysis/analysis.md produced by `run.py inspect`).
"""
import numpy as np
from escher_handlayout.engine import unit


def build(L):
    # ---- 1. subsystems: name, line colour, title shown on the map, text of the EXISTING title label to reuse (or None)
    L.zone('Subsystem A', '#2b5c9e', 'Subsystem A', 'subsystem a label text in the input map')
    # L.zone('Subsystem B', '#5b8c2a', 'Subsystem B', None)

    # titles of pieces that may end up disconnected (optional; defaults are generated automatically)
    # L.small_titles['SOME_RXN_BIGG_ID'] = ('Subsystem A', 'short title', None)

    L.bend = 45          # shift cofactor dots along the axis so that their curves have room (45..70)

    # ---- 2. splits and merges
    # copy = L.split('<primary node id>', <reaction id>)   # a second node of that metabolite used only by that reaction
    # L.merge('<keep id>', '<dup id 1>', '<dup id 2>')     # one hub node for the same metabolite in several subsystems

    # ---- 3. coordinates of primary metabolites + 4. reactions of subsystem A
    # L.setp({101: (0, 0), 102: (0, 300), 103: (0, 600)})
    # L.reaction(11, 'Subsystem A', side='R')                       # vertical chain step, cofactors on the right
    # L.reaction(12, 'Subsystem A', shift=(-110, 0), side='L')      # reverse reaction as a lane on the left
    # L.reaction(13, 'Subsystem A', axis=(0, 1), mid=(0, 450), lab='L')   # branch point: vertical axis, custom midpoint
    # ring (cycle):  L.reaction(rid, zone, arc=((cx, cy), radial_offset), side='out'|'in')
    # parallel isozymes on one lane: L.reaction(rid, zone, shift=(0, dy), side='R', curve=False, gap=40)

    # ---- 5. title position of every subsystem that keeps at least one reaction in the main component
    # L.titles['Subsystem A'] = (x, y)
