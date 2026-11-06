"""Drawing the maps: a box per subsystem with two zoom levels, the links between the boxes, the view; maps with an
AI layout as their Escher drawing with the outlines of the subsystems."""

from __future__ import annotations

import re

import numpy as np
from qtpy.QtCore import QByteArray, QLineF, QPointF, QRectF, Qt, Signal
from qtpy.QtGui import QColor, QFont, QFontMetricsF, QPainter, QPen
from qtpy.QtSvg import QSvgRenderer
from qtpy.QtWidgets import (
    QGraphicsItem,
    QGraphicsLineItem,
    QGraphicsObject,
    QGraphicsScene,
    QGraphicsView,
    QStyleOptionGraphicsItem,
)

from .session import Session, boxes_free, constrained_delta

#: zoom (screen pixels per map unit) from which the subsystems are drawn as SVG in the mode AUTO
DETAIL_SCALE = 0.3
AUTO, BOXES, DRAWINGS = "auto", "boxes", "drawings"
SELECTED = QColor(30, 110, 230)
FOUND = QColor(240, 160, 0)
#: reactions added from the model (not on the BRENDA map): the stroke of their enzyme nodes in the SVG
ADDED = QColor(130, 60, 160)
TEXT = QColor(25, 25, 25)


def qcolor(text: str) -> QColor:
    """A colour of the BRENDA map ("rgb(r, g, b)" or a name)."""
    m = re.match(r"\s*rgba?\(\s*(\d+)\s*,\s*(\d+)\s*,\s*(\d+)", text or "")
    return QColor(*(int(v) for v in m.groups())) if m else QColor(text or "#dedede")


class SubsystemItem(QGraphicsObject):
    """A subsystem on a map: its box with name (zoomed out) or its SVG drawing (zoomed in); on a map shown as its
    AI layout, the outline of its nodes there (`outline`: x, y, w, h).

    The local coordinates are those of the box (top left at 0, 0, map
    units); the title of the drawing sticks out above the box.
    """

    def __init__(self, map_scene: MapScene, key: str, outline=None):
        super().__init__()
        self.map_scene, self.key = map_scene, key
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        self.found = False
        self.renderer = None
        self.outline = outline
        self.refresh()

    def refresh(self) -> None:
        """Size, texts and drawing from the session (after changes of the assignment)."""
        session = self.map_scene.session
        s = session.dm.subsystems[self.key]
        self.prepareGeometryChange()
        self.w, self.h, self.label, self.fill = s.w, s.h, s.label, qcolor(s.fill)
        self.renderer = None
        if self.outline is not None:
            self.w, self.h = self.outline[2], self.outline[3]
            self.view_box = self.box()
        else:
            self.view_box = QRectF(*session.subsystem_svg(self.map_scene.name, self.key)[1])
        others = [m for m in session.maps_of(self.key) if m != self.map_scene.name]
        count, added = session.counts.get(self.key, 0), session.added_counts.get(self.key, 0)
        self.lines = [(self.label, True, 15), (f"{count} reactions", False, 12)]
        if added:
            self.lines.append((f"{added} added from the model", False, 12, ADDED))
        if others:
            self.lines.append(("also on: " + ", ".join(others), False, 11))
        self.setToolTip(f"{self.label} ({self.key})\n{count} reactions"
                        + (f", {added} added from the model" if added else "")
                        + (f"\nalso on: {', '.join(others)}" if others else ""))
        self.update()

    def box(self) -> QRectF:
        return QRectF(0, 0, self.w, self.h)

    def boundingRect(self) -> QRectF:
        return self.view_box.united(self.box().adjusted(-4, -4, 4, 4))

    def detailed(self, scale: float) -> bool:
        mode = self.map_scene.mode
        return mode == DRAWINGS or (mode == AUTO and scale >= DETAIL_SCALE)

    def paint(self, painter: QPainter, option, widget=None) -> None:
        scale = QStyleOptionGraphicsItem.levelOfDetailFromTransform(painter.worldTransform())
        box = self.box()
        if self.outline is not None:
            marked = self.isSelected() or self.found
            pen = QPen(SELECTED if self.isSelected() else FOUND if self.found else self.fill.darker(170),
                       4 if marked else 1.5, Qt.PenStyle.SolidLine if marked else Qt.PenStyle.DashLine)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(box)
            return
        if self.detailed(scale):
            if self.renderer is None:
                data, _ = self.map_scene.session.subsystem_svg(self.map_scene.name, self.key)
                self.renderer = QSvgRenderer(QByteArray(data))
            self.renderer.render(painter, self.view_box)
        else:
            fill = QColor(self.fill)
            fill.setAlphaF(0.6)
            pen = QPen(self.fill.darker(150), 1.5)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.setBrush(fill)
            painter.drawRoundedRect(box, 18, 18)
            self.paint_texts(painter, box, scale)
        if self.isSelected() or self.found:
            pen = QPen(SELECTED if self.isSelected() else FOUND, 4 if self.isSelected() else 3)
            pen.setCosmetic(True)
            painter.setPen(pen)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(box)

    def paint_texts(self, painter: QPainter, box: QRectF, scale: float) -> None:
        """Name, number of reactions and other maps, in a fixed size on the screen (cut to the box)."""
        width, height = box.width() * scale - 8, box.height() * scale - 4
        if width < 24 or height < 10:
            return
        fonts = []
        for text, bold, size, *colour in self.lines:
            font = QFont()
            font.setPixelSize(size)
            font.setBold(bold)
            fonts.append((text, font, colour[0] if colour else TEXT, QFontMetricsF(font)))
        total = sum(fm.height() for *_, fm in fonts)
        painter.save()
        painter.setClipRect(box)
        painter.translate(box.center())
        painter.scale(1 / scale, 1 / scale)
        y = -min(total, height) / 2
        for text, font, colour, fm in fonts:
            if y + fm.height() > height / 2 + 1:
                break
            painter.setFont(font)
            painter.setPen(colour)
            painter.drawText(QRectF(-width / 2, y, width, fm.height()), Qt.AlignmentFlag.AlignCenter,
                             fm.elidedText(text, Qt.TextElideMode.ElideRight, width))
            y += fm.height()
        painter.restore()


class EscherItem(QGraphicsObject):
    """The AI layout of a map: the SVG export of its Escher map in map units (or a note if there is none)."""

    def __init__(self, layout):
        super().__init__()
        self.setZValue(-1)
        self.setCacheMode(QGraphicsItem.CacheMode.DeviceCoordinateCache)
        self.set_layout(layout)

    def set_layout(self, layout) -> None:
        self.prepareGeometryChange()
        self.layout = layout
        self.rect = QRectF(*layout.extent())
        self.renderer = QSvgRenderer(QByteArray(layout.svg.encode("utf-8"))) if layout.svg else None
        self.update()

    def boundingRect(self) -> QRectF:
        return self.rect

    def paint(self, painter: QPainter, option, widget=None) -> None:
        if self.renderer is not None and self.renderer.isValid():
            self.renderer.render(painter, self.rect)
            return
        painter.fillRect(self.rect, QColor(255, 255, 255))
        scale = QStyleOptionGraphicsItem.levelOfDetailFromTransform(painter.worldTransform())
        painter.save()
        painter.translate(self.rect.center())
        painter.scale(1 / scale, 1 / scale)
        font = QFont()
        font.setPixelSize(14)
        painter.setFont(font)
        painter.setPen(TEXT)
        painter.drawText(QRectF(-250, -30, 500, 60), Qt.AlignmentFlag.AlignCenter,
                         "AI layout without SVG export\n(Escher could not be run: no QtWebEngine or escher.min.js)")
        painter.restore()


class MapScene(QGraphicsScene):
    """The subsystems of one map (`name`; None: the subsystems on no map, which can not be moved).

    A map with an AI layout is shown as that (`ai_view`, :class:`EscherItem`)
    with the outlines of its subsystems, which can be selected, not moved.
    """

    def __init__(self, session: Session, name: str | None, parent=None):
        super().__init__(parent)
        self.session, self.name = session, name
        self.boxes: dict[str, SubsystemItem] = {}
        self.link_items: list[tuple[QGraphicsLineItem, str, QPointF, str, QPointF]] = []
        self.mode = AUTO
        self.version = None
        self.ai_view = True
        self.escher_item = None
        self.shown = None  # the layout shown (by the items)
        self.setBackgroundBrush(QColor(250, 250, 250))

    def ai_info(self):
        """(layout, outlines of the subsystems, problems) if the map is shown as its AI layout, else None."""
        return self.session.layout_info(self.name) if self.ai_view and self.name is not None else None

    @property
    def read_only(self) -> bool:
        return self.name is None or self.ai_info() is not None

    def positions(self) -> dict[str, tuple[float, float]]:
        if self.name is None:
            return self.session.unassigned_positions()
        info = self.ai_info()
        if info is not None:
            return {k: (x, y) for k, (x, y, _, _) in info[1].items() if k in self.session.positions[self.name]}
        return self.session.positions[self.name]

    def sync(self) -> None:
        """Items and positions as in the session; the drawings again if the assignment changed."""
        positions = self.positions()
        info = self.ai_info()
        layout = info[0] if info is not None else None
        redraw = self.version != self.session.version or self.shown is not layout
        if layout is None and self.escher_item is not None:
            self.removeItem(self.escher_item)
            self.escher_item = None
        elif layout is not None and self.escher_item is None:
            self.escher_item = EscherItem(layout)
            self.addItem(self.escher_item)
        elif layout is not None and self.escher_item.layout is not layout:
            self.escher_item.set_layout(layout)
        for k in [k for k in self.boxes if k not in positions]:
            self.removeItem(self.boxes.pop(k))
        for k, (x, y) in positions.items():
            outline = info[1][k] if info is not None else None
            item = self.boxes.get(k)
            if item is None:
                item = self.boxes[k] = SubsystemItem(self, k, outline)
                self.addItem(item)
            elif redraw:
                item.outline = outline
                item.refresh()
            item.setPos(x, y)
        self.version, self.shown = self.session.version, layout
        if redraw:
            self.make_links()
        self.update_links()
        self.fit_scene_rect()

    def set_ai_view(self, on: bool) -> None:
        self.ai_view = on
        self.sync()

    def make_links(self) -> None:
        for line, *_ in self.link_items:
            self.removeItem(line)
        self.link_items = []
        if self.name is None or self.ai_info() is not None:
            return
        pen = QPen(QColor(120, 120, 120, 170), 2.5, Qt.PenStyle.DashLine)
        pen.setCosmetic(True)
        for a, pa, b, pb in self.session.links(self.name):
            line = QGraphicsLineItem()
            line.setPen(pen)
            line.setZValue(1)
            self.addItem(line)
            self.link_items.append((line, a, QPointF(*pa), b, QPointF(*pb)))

    def update_links(self) -> None:
        for line, a, pa, b, pb in self.link_items:
            if a in self.boxes and b in self.boxes:
                line.setLine(QLineF(self.boxes[a].pos() + pa, self.boxes[b].pos() + pb))

    def fit_scene_rect(self) -> None:
        """Room around the boxes to move them into."""
        rect = self.itemsBoundingRect()
        if rect.isEmpty():
            rect = QRectF(0, 0, 1000, 1000)
        m = max(rect.width(), rect.height(), 2000.0)
        self.setSceneRect(rect.adjusted(-m, -m, m, m))

    def content_rect(self) -> QRectF:
        rect = QRectF()
        for item in [*self.boxes.values(), *([self.escher_item] if self.escher_item is not None else [])]:
            rect = rect.united(item.mapRectToScene(item.boundingRect()))
        return rect

    def selected_keys(self) -> list[str]:
        return [item.key for item in self.selectedItems() if isinstance(item, SubsystemItem)]

    def select(self, keys) -> None:
        self.clearSelection()
        for k in keys:
            if k in self.boxes:
                self.boxes[k].setSelected(True)

    def set_found(self, keys) -> int:
        keys = set(keys)
        for k, item in self.boxes.items():
            if item.found != (k in keys):
                item.found = k in keys
                item.update()
        return len(keys & set(self.boxes))

    def set_mode(self, mode: str) -> None:
        self.mode = mode
        self.update()


class MapView(QGraphicsView):
    """A map: zoom with the wheel, pan with the middle button, select (click, Ctrl+click, rubber band) and drag.

    The selected boxes are dragged together and stop at (or slide along)
    other boxes, at least the gap from them (:func:`.session.constrained_delta`).
    Signals: `moved({subsystem ID: (x, y)})` after a drag, `menu(global
    position)` for the context menu, `zoom_changed()`.
    """

    moved = Signal(object)
    menu = Signal(object)
    zoom_changed = Signal()

    def __init__(self, scene: MapScene, parent=None):
        super().__init__(scene, parent)
        self.setRenderHints(QPainter.RenderHint.Antialiasing | QPainter.RenderHint.TextAntialiasing
                            | QPainter.RenderHint.SmoothPixmapTransform)
        self.setDragMode(QGraphicsView.DragMode.RubberBandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setViewportUpdateMode(QGraphicsView.ViewportUpdateMode.SmartViewportUpdate)
        self._drag = None
        self._pan = None
        self.fitted = False

    def map_scene(self) -> MapScene:
        return self.scene()

    def scale_factor(self) -> float:
        return self.transform().m11()

    def zoom(self, factor: float) -> None:
        factor = min(max(factor, 0.002 / self.scale_factor()), 6.0 / self.scale_factor())
        self.scale(factor, factor)
        self.zoom_changed.emit()

    def zoom_to(self, rect: QRectF) -> None:
        if not rect.isEmpty():
            m = 0.04 * max(rect.width(), rect.height())
            self.fitInView(rect.adjusted(-m, -m, m, m), Qt.AspectRatioMode.KeepAspectRatio)
            self.zoom_changed.emit()

    def zoom_to_fit(self) -> None:
        self.zoom_to(self.map_scene().content_rect())
        self.fitted = True

    def showEvent(self, event) -> None:
        super().showEvent(event)
        if not self.fitted:
            self.zoom_to_fit()

    def subsystem_at(self, pos) -> SubsystemItem | None:
        for item in self.items(pos):
            if isinstance(item, SubsystemItem):
                return item
        return None

    def wheelEvent(self, event) -> None:
        self.zoom(1.0015 ** event.angleDelta().y())

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.MiddleButton:
            self._pan = event.position()
            self.viewport().setCursor(Qt.CursorShape.ClosedHandCursor)
            return
        item = self.subsystem_at(event.position().toPoint()) if event.button() == Qt.MouseButton.LeftButton else None
        if item is not None:
            scene = self.map_scene()
            if event.modifiers() & Qt.KeyboardModifier.ControlModifier:
                item.setSelected(not item.isSelected())
                return
            if not item.isSelected():
                scene.clearSelection()
                item.setSelected(True)
            if not scene.read_only:
                session, keys = scene.session, scene.selected_keys()
                start = self.mapToScene(event.position().toPoint())
                self._drag = dict(start=start, keys=keys, origin={k: scene.boxes[k].pos() for k in keys},
                                  moving=np.array([session.rect(scene.name, k) for k in keys]),
                                  others=session.rects(scene.name, exclude=keys), offset=(0.0, 0.0))
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:
        if self._pan is not None:
            delta = event.position() - self._pan
            self._pan = event.position()
            self.horizontalScrollBar().setValue(int(self.horizontalScrollBar().value() - delta.x()))
            self.verticalScrollBar().setValue(int(self.verticalScrollBar().value() - delta.y()))
            return
        if self._drag is not None:
            d, scene = self._drag, self.map_scene()
            p = self.mapToScene(event.position().toPoint())
            tx, ty = p.x() - d["start"].x(), p.y() - d["start"].y()
            gap = scene.session.gap
            if boxes_free(d["moving"] + [tx, ty, tx, ty], d["others"], gap):
                offset = (tx, ty)  # free at the cursor (also behind an obstacle)
            else:
                ox, oy = d["offset"]
                ddx, ddy = constrained_delta(d["moving"] + [ox, oy, ox, oy], d["others"], tx - ox, ty - oy, gap)
                offset = (ox + ddx, oy + ddy)
            d["offset"] = offset
            for k in d["keys"]:
                scene.boxes[k].setPos(d["origin"][k] + QPointF(*offset))
            scene.update_links()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:
        if self._pan is not None and event.button() == Qt.MouseButton.MiddleButton:
            self._pan = None
            self.viewport().unsetCursor()
            return
        if self._drag is not None and event.button() == Qt.MouseButton.LeftButton:
            d, self._drag = self._drag, None
            ox, oy = d["offset"]
            if abs(ox) > 0.01 or abs(oy) > 0.01:
                self.moved.emit({k: (p.x() + ox, p.y() + oy) for k, p in d["origin"].items()})
            return
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:
        item = self.subsystem_at(event.position().toPoint())
        if item is not None:
            self.zoom_to(item.mapRectToScene(item.boundingRect()))
        else:
            super().mouseDoubleClickEvent(event)

    def contextMenuEvent(self, event) -> None:
        item = self.subsystem_at(event.pos())
        if item is not None and not item.isSelected():
            self.map_scene().clearSelection()
            item.setSelected(True)
        self.menu.emit(event.globalPos())

    def cursor_scene_pos(self, global_pos) -> tuple[float, float] | None:
        """The scene position under the mouse pointer if it is over the view."""
        p = self.viewport().mapFromGlobal(global_pos)
        if not self.viewport().rect().contains(p):
            return None
        s = self.mapToScene(p)
        return s.x(), s.y()
