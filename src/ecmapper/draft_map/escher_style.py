"""Escher maps of draft maps in the style of the Escher maps of the BRENDA map (:mod:`ecmapper.brenda_escher`).

The draft style (:func:`.escher.to_escher`) keeps every node where the draft
map has it. The BRENDA style draws each subsystem as :func:`ecmapper.brenda_escher.build`
does (rows and columns, rounded corners, cofactors on branches), in its place
on the map. Both record the draft objects of the Escher IDs alike, so AI
layouts (:mod:`.ai_layouts`) work with either.
"""

from __future__ import annotations

import re
from collections import defaultdict

from .core import DraftMap, ModelData
from .escher import SCHEMA, _Ids, place_text_labels, stub_text
from .relayout import subsystem_map

#: the styles of exported Escher maps: {style: description}
STYLES = {"draft": "as the draft map (nodes where the draft map has them)",
          "brenda": "as the Escher maps of the BRENDA map (ecmapper.brenda_escher)"}


def subsystem_escher(sub: DraftMap, md: ModelData, label: str, scale: float = 2.0, new_id=None) -> list | None:
    """The Escher map of a draft map of one subsystem (:func:`.relayout.subsystem_map`) in the BRENDA style, with its
    title; None if it has no model reactions."""
    from ..brenda_escher import build, from_draft

    reactions, data = from_draft(sub, md, scale)
    if not reactions:
        return None
    return build(reactions, data, label, description="", new_id=new_id)


def to_brenda_escher(m: DraftMap, md: ModelData, map_name: str, scale: float = 2.0, subsystem_labels: bool = True,
                     description: str = "", ids: dict | None = None, reuse: dict | None = None) -> list:
    """The Escher map of map `m` (absolute positions) in the BRENDA style: each subsystem drawn by
    :func:`subsystem_escher` in its place, the "→ map" labels of metabolites that continue on other maps next to
    them. `ids` and `reuse`: as :func:`.escher.to_escher`."""
    ids = {} if ids is None else ids  # the stubs are found by it
    new_id = _Ids(reuse, ids)
    nodes, reactions, text_labels = {}, {}, {}
    for k, s in m.subsystems.items():
        sub = subsystem_map(m, k)
        data = subsystem_escher(sub, md, s.label or k, scale, new_id)
        if data is None:
            continue
        body = data[1]
        nodes.update(body["nodes"])
        reactions.update(body["reactions"])
        if subsystem_labels:
            text_labels.update(body["text_labels"])
    at = defaultdict(list)  # drawn node -> its Escher nodes
    for i, key in ids.items():
        if key[0] == "node" and i in nodes:
            at[key[1]].append(nodes[i])
    references = {}
    for drawn_id, escher_nodes in at.items():
        drawn = m.nodes.get(drawn_id)
        if drawn is None or not getattr(drawn, "stubs", None):
            continue
        xs, ys = [n["x"] for n in escher_nodes], [n["y"] for n in escher_nodes]
        label_id = new_id(("stub", drawn_id))
        text_labels[label_id] = dict(text=stub_text(drawn.stubs), x=xs[0], y=max(ys) + 66)
        references[label_id] = (sum(xs) / len(xs), min(ys), max(ys))
    body = dict(reactions=reactions, nodes=nodes, text_labels=text_labels)
    if references:
        place_text_labels(body, references)
    xs = [n["x"] for n in nodes.values()] + [t["x"] for t in text_labels.values()]
    xs += [t["x"] + 28.0 * len(t["text"]) for t in text_labels.values()]
    ys = [n["y"] for n in nodes.values()] + [t["y"] for t in text_labels.values()]
    if not xs:
        xs, ys = [0.0], [0.0]
    margin = 250
    body["canvas"] = dict(x=min(xs) - margin, y=min(ys) - margin, width=max(xs) - min(xs) + 2 * margin,
                          height=max(ys) - min(ys) + 2 * margin)
    head = dict(map_name=map_name, map_id=re.sub(r"\W+", "_", map_name), map_description=description,
                homepage="https://escher.github.io", schema=SCHEMA)
    return [head, body]
