"""Lexer, parser and document model of PWL 2.0 (see ``PWL.md`` and ``pwl.ebnf``).

The parser has no dependencies. :func:`parse` returns a :class:`Document`;
syntax errors raise :class:`PWLError` with line and column. :meth:`Document.problems`
lists the semantic errors (reaction IDs used twice, unknown reaction IDs in
motifs, ...), :func:`parse` raises them as well unless ``validate=False``.

PWL is primarily a layout description: compound declarations, reaction IDs,
enzymes, coefficients and directions are optional. ``direction`` statements
set or change directions after the layout statements.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterator

VERSION = "2.0"

#: words that start a statement or end a chain; compound tokens must not be one of them
KEYWORDS = frozenset(
    {"pwl", "subsystem", "compartment", "compound", "flow", "port", "ring", "chord", "star", "react", "motif",
     "effector", "wrap", "direction"}
)
FLOWS = ("down", "up", "left", "right")
PORT_SIDES = ("U", "D", "L", "R")
SENSES = ("cw", "ccw")
WRAPS = ("serpentine", "spiral", "column")
EFFECTS = ("inhibits", "activates")
#: link arrows (opening part, closing part) -> direction of the reaction from the left to the right node
LINK_ARROWS = {("-", "-"): None, ("-", "->"): "->", ("<-", "-"): "<-", ("<-", "->"): "<->", ("-", "=>"): "=>",
               ("<=", "-"): "<="}
DIRECTIONS = ("->", "<-", "<->", "=>", "<=")
#: motif kinds with the type of their members ("compound" tokens or "reaction" IDs) and their roles
MOTIF_KINDS = {
    "cycle": ("compound", {}),
    "substrate_cycle": ("compound", {"forward": "reaction", "backward": "reaction"}),
    "shunt": ("compound", {"bypasses": "compound"}),
    "repeat": ("reaction", {}),
    "carrier": ("compound", {"forward": "reaction", "backward": "reaction"}),
    "transport": ("compound", {"via": "reaction"}),
}


class PWLError(ValueError):
    """Syntax error in a PWL text."""

    def __init__(self, message: str, line: int = 0, column: int = 0):
        super().__init__(f"line {line}:{column}: {message}" if line else message)
        self.line, self.column = line, column


# ---------------------------------------------------------------------------------- model


@dataclass
class Participant:
    """A compound in a reaction: ``[coef*]token[#copy][@compartment]``.

    `copy` numbers the drawn copies of a compound in a subsystem (``pyruvate#2``).
    """

    token: str
    copy: int = 1
    coef: str | None = None
    location: str | None = None

    @property
    def node(self) -> str:
        """The drawn node: token and copy number."""
        return self.token if self.copy == 1 else f"{self.token}#{self.copy}"

    @property
    def species(self) -> str:
        """The compound in its compartment, e.g. ``H@side2``."""
        return self.token if self.location is None else f"{self.token}@{self.location}"

    def __str__(self) -> str:
        text = self.node if self.location is None else f"{self.node}@{self.location}"
        return text if self.coef is None else f"{self.coef}*{text}"


@dataclass
class Node:
    """A position in a chain: one participant, a group ``(A + B)`` of the same side, or the open end ``_``."""

    participants: list[Participant]
    group: bool = False

    @property
    def open(self) -> bool:
        return not self.participants

    def __str__(self) -> str:
        if self.open:
            return "_"
        if self.group:
            return "(" + " + ".join(map(str, self.participants)) + ")"
        return str(self.participants[0])


@dataclass
class CofactorItem:
    """``a>b`` (couple: a consumed, b produced), ``-a`` (consumed) or ``+a`` (produced)."""

    kind: str  # couple, consumed, produced
    a: Participant
    b: Participant | None = None

    def __str__(self) -> str:
        if self.kind == "couple":
            return f"{self.a}>{self.b}"
        return ("-" if self.kind == "consumed" else "+") + str(self.a)


@dataclass
class Reaction:
    """One reaction (enzyme node).

    `left` and `right` are the participants written before and after the link
    (main metabolites); for a ``react`` statement they are empty and the
    participants are in `unsided`. `direction` is one of :data:`DIRECTIONS`
    (from left to right) or None if the direction is not known.
    """

    enzymes: list[str]
    spontaneous: bool = False
    rid: str | None = None
    cofactors: list[CofactorItem] = field(default_factory=list)
    cofactor_side: str | None = None
    attrs: dict[str, str] = field(default_factory=dict)
    direction: str | None = None
    left: list[Participant] = field(default_factory=list)
    right: list[Participant] = field(default_factory=list)
    unsided: list[Participant] = field(default_factory=list)
    line: int = 0
    comment: str = ""
    context: str = "chain"  # chain, ring, chord, star, react

    @property
    def consumed(self) -> list[Participant]:
        return [c.a for c in self.cofactors if c.kind in ("couple", "consumed")]

    @property
    def produced(self) -> list[Participant]:
        return [c.b if c.kind == "couple" else c.a for c in self.cofactors if c.kind in ("couple", "produced")]

    @property
    def sided(self) -> bool:
        """Whether the main metabolites are split into two sides of the reaction."""
        return not self.unsided

    @property
    def substrates(self) -> list[Participant] | None:
        """Main substrates in the direction of the reaction, None if the direction is not known."""
        if self.direction in ("->", "=>", "<->"):
            return self.left
        if self.direction in ("<-", "<="):
            return self.right
        return None

    @property
    def products(self) -> list[Participant] | None:
        if self.direction in ("->", "=>", "<->"):
            return self.right
        if self.direction in ("<-", "<="):
            return self.left
        return None

    @property
    def reversible(self) -> bool | None:
        """True for ``<->``, False for ``=>``/``<=``, None if not stated."""
        if self.direction == "<->":
            return True
        if self.direction in ("=>", "<="):
            return False
        return None

    @property
    def participants(self) -> list[Participant]:
        return self.left + self.right + self.unsided

    @property
    def enzyme_label(self) -> str:
        """EC numbers joined by '/', '?' for an unknown enzyme, '' for a reaction given only by its ID."""
        return "/".join(self.enzymes)


@dataclass
class Link:
    reactions: list[Reaction]
    turn: str = ""
    arrow: str | None = None  # direction written on the link, None for -[..]-
    line: int = 0


@dataclass
class Chain:
    nodes: list[Node]
    links: list[Link]
    chord: bool = False
    wrap: tuple[str, int] | None = None
    line: int = 0


@dataclass
class Ring:
    sense: str
    statements: list  # Chain and React
    line: int = 0


@dataclass
class Star:
    hub: Node
    arc: int
    order: str
    spokes: list[Chain]  # each spoke starts at the hub
    line: int = 0


@dataclass
class React:
    reaction: Reaction
    line: int = 0


@dataclass
class Compound:
    token: str
    label: str | None = None
    xrefs: list[str] = field(default_factory=list)
    line: int = 0

    def xref(self, namespace: str) -> list[str]:
        """The identifiers of the compound in a namespace, e.g. ``brenda``."""
        return [x.split(":", 1)[1] for x in self.xrefs if x.split(":", 1)[0] == namespace]


@dataclass
class Port:
    participant: Participant
    side: str
    target: str | None = None
    line: int = 0


@dataclass
class Motif:
    kind: str
    members: list[str]
    params: dict[str, str] = field(default_factory=dict)
    roles: dict[str, list[str]] = field(default_factory=dict)
    line: int = 0


@dataclass
class Effector:
    participant: Participant
    effect: str
    reactions: list[str]
    line: int = 0


@dataclass
class DirectionSetting:
    """``direction r1, r2 ->;``: sets (or changes) the direction of reactions after the layout statements.

    The direction is relative to the order in which the reaction's sides are
    written (left node to right node); ``?`` makes it unknown again.
    """

    reactions: list[str]
    direction: str | None
    line: int = 0


@dataclass
class Document:
    id: str
    title: str
    version: str = VERSION
    flow: str | None = None
    compartments: dict[str, str | None] = field(default_factory=dict)
    compounds: dict[str, Compound] = field(default_factory=dict)
    ports: list[Port] = field(default_factory=list)
    statements: list = field(default_factory=list)  # Chain, Ring, Star, React
    motifs: list[Motif] = field(default_factory=list)
    effectors: list[Effector] = field(default_factory=list)
    directions: list[DirectionSetting] = field(default_factory=list)

    def reactions(self) -> list[Reaction]:
        """All reactions in the order of the text."""
        return list(_walk_reactions(self.statements))

    def reaction_index(self) -> dict[str, Reaction]:
        return {r.rid: r for r in self.reactions() if r.rid}

    def nodes(self) -> set[str]:
        """The drawn nodes (``token#copy``) used by the reactions."""
        return {p.node for r in self.reactions() for p in r.participants}

    def problems(self) -> list[str]:
        """Semantic errors, see :func:`validate`."""
        return validate(self)


def _walk_reactions(statements) -> Iterator[Reaction]:
    for s in statements:
        if isinstance(s, Ring):
            yield from _walk_reactions(s.statements)
        elif isinstance(s, Star):
            for spoke in s.spokes:
                for link in spoke.links:
                    yield from link.reactions
        elif isinstance(s, Chain):
            for link in s.links:
                yield from link.reactions
        elif isinstance(s, React):
            yield s.reaction


# ---------------------------------------------------------------------------------- lexer

_STR = re.compile(r"'(?:[^'\\\n]|\\.)*'")
_BRACKET = re.compile(r"(?:(<=|<-|-)([FLRB])?)?\[([^\]\n]*)\](=>|->|-)?")
_NUM = re.compile(r"\d+(?:\.\d+)?(?![A-Za-z0-9_.])")
_ID = re.compile(r"[A-Za-z0-9_]+")
_WIDE_ID = re.compile(r"[A-Za-z0-9_][A-Za-z0-9_\-]*")
_PUNCT = set("{}();,+*#@:=?")
_DIR = re.compile(r"<->|->|<-|=>|<=")


@dataclass
class _Tok:
    kind: str  # ID, NUM, STR, BRK, PUNCT, EOF
    value: object
    line: int
    col: int


def _unquote(s: str) -> str:
    return re.sub(r"\\(.)", r"\1", s[1:-1])


def quote(s: str) -> str:
    """A PWL string literal."""
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def _lex(text: str) -> tuple[list[_Tok], dict[int, str]]:
    toks: list[_Tok] = []
    comments: dict[int, str] = {}
    for ln, line in enumerate(text.split("\n"), 1):
        pos = 0
        n = len(line)
        while pos < n:
            c = line[pos]
            if c.isspace():
                pos += 1
                continue
            if line.startswith("//", pos):
                comments[ln] = line[pos + 2:].strip()
                break
            col = pos + 1
            if c == "'":
                m = _STR.match(line, pos)
                if not m:
                    raise PWLError("unterminated string", ln, col)
                toks.append(_Tok("STR", _unquote(m.group(0)), ln, col))
                pos = m.end()
                continue
            m = _BRACKET.match(line, pos)
            if m and (c in "-<[" ):
                toks.append(_Tok("BRK", (m.group(1), m.group(2) or "", m.group(3), m.group(4)), ln, col))
                pos = m.end()
                continue
            m = _DIR.match(line, pos)
            if m and not line.startswith("[", m.end()):
                toks.append(_Tok("DIR", m.group(0), ln, col))
                pos = m.end()
                continue
            if c in "-<[":
                raise PWLError("malformed link, expected -[...]-", ln, col)
            # subsystem IDs may contain '-': only directly after 'subsystem' and 'to'
            if toks and toks[-1].kind == "ID" and toks[-1].value in ("subsystem", "to"):
                m = _WIDE_ID.match(line, pos)
                if m:
                    toks.append(_Tok("ID", m.group(0), ln, col))
                    pos = m.end()
                    continue
            m = _NUM.match(line, pos)
            if m:
                toks.append(_Tok("NUM", m.group(0), ln, col))
                pos = m.end()
                continue
            m = _ID.match(line, pos)
            if m:
                toks.append(_Tok("ID", m.group(0), ln, col))
                pos = m.end()
                continue
            if c in _PUNCT:
                toks.append(_Tok("PUNCT", c, ln, col))
                pos += 1
                continue
            raise PWLError(f"unexpected character {c!r}", ln, col)
    toks.append(_Tok("EOF", "", len(text.split("\n")), 0))
    return toks, comments


# lexer of the reaction bodies inside [...]
_BODY = re.compile(
    r"\s*(?:(?P<dir><->|->|<-|=>|<=)"
    r"|(?P<ec>\d+(?:\.[0-9A-Za-z?\-]+){2,3}(?:/\d+(?:\.[0-9A-Za-z?\-]+){2,3})*)"
    r"|(?P<num>\d+(?:\.\d+)?(?![A-Za-z0-9_.]))"
    r"|(?P<id>[A-Za-z0-9_]+)"
    r"|(?P<punct>[|,>+\-*@^:~?]))"
)
_ATTRS = re.compile(r"\{([^{}]*)\}\s*$")
_ATTR = re.compile(r"\s*([A-Za-z_][A-Za-z0-9_]*)\s*=\s*([A-Za-z0-9_.+\-:/]+)\s*$")


def _lex_body(text: str, line: int, col: int) -> list[tuple[str, str]]:
    out = []
    pos = 0
    text = text.rstrip()
    while pos < len(text):
        m = _BODY.match(text, pos)
        if not m or m.end() == pos:
            raise PWLError(f"unexpected {text[pos:].strip()[:10]!r} in reaction", line, col + pos)
        kind = m.lastgroup
        out.append((kind, m.group(kind)))
        pos = m.end()
    return out


def _split_reactions(body: str) -> list[str]:
    return [part for part in body.split("||")]


def _parse_reaction(text: str, line: int, col: int) -> Reaction:
    attrs: dict[str, str] = {}
    m = _ATTRS.search(text)
    if m:
        for item in m.group(1).split(","):
            if not item.strip():
                continue
            a = _ATTR.match(item)
            if not a:
                raise PWLError(f"bad attribute {item.strip()!r}", line, col)
            if a.group(1) in attrs:
                raise PWLError(f"attribute {a.group(1)!r} given twice", line, col)
            attrs[a.group(1)] = a.group(2)
        text = text[: m.start()]
    toks = _lex_body(text, line, col)
    i = 0

    def peek(k=0):
        return toks[i + k] if i + k < len(toks) else ("end", "")

    def take(kind=None, value=None):
        nonlocal i
        t = peek()
        if (kind and t[0] != kind) or (value is not None and t[1] != value):
            want = value or kind
            raise PWLError(f"expected {want} in reaction, got {t[1] or 'end'!r}", line, col)
        i += 1
        return t

    rid = None
    if peek()[0] == "id" and peek(1) == ("punct", ":"):
        rid = take()[1]
        take()
    direction = take()[1] if peek()[0] == "dir" else None
    spont = False
    if peek() == ("punct", "~"):
        take()
        spont = True
    if peek() == ("punct", "?"):
        take()
        enzymes = ["?"]
    elif peek()[0] == "ec":
        enzymes = take()[1].split("/")
    elif peek()[0] == "id" and rid is None and not spont and peek(1) not in (("punct", "*"), ("punct", ">")):
        # a reaction given only by its ID (any name, no primary data source needed): [PGI | ATP>ADP]
        rid = take()[1]
        enzymes = []
    elif peek()[0] in ("end",) or peek() in (("punct", "|"), ("punct", "^")):
        enzymes = []  # no enzyme: [r1:], [], [| ATP>ADP]
    else:
        raise PWLError(f"expected EC number, '?' or a reaction ID in reaction, got {peek()[1]!r}", line, col)
    cofactors = []

    def participant() -> Participant:
        coef = None
        if peek()[0] in ("num", "id") and peek(1) == ("punct", "*"):
            coef = take()[1]
            take()
            if not (_NUM.fullmatch(coef) or coef == "n"):
                raise PWLError(f"bad coefficient {coef!r}", line, col)
        token = take("id")[1]
        loc = None
        if peek() == ("punct", "@"):
            take()
            loc = take("id")[1]
        return Participant(token, 1, coef, loc)

    if peek() == ("punct", "|"):
        take()
        while True:
            if peek() in (("punct", "-"), ("punct", "+")):
                sign = take()[1]
                cofactors.append(CofactorItem("consumed" if sign == "-" else "produced", participant()))
            else:
                a = participant()
                take("punct", ">")
                cofactors.append(CofactorItem("couple", a, participant()))
            if peek() == ("punct", ","):
                take()
                continue
            break
    side = None
    if peek() == ("punct", "^"):
        take()
        side = take("id")[1]
        if side not in ("L", "R"):
            raise PWLError(f"cofactor side must be L or R, not {side!r}", line, col)
        if not cofactors:
            raise PWLError("cofactor side without cofactors", line, col)
    if i != len(toks):
        raise PWLError(f"unexpected {peek()[1]!r} in reaction", line, col)
    return Reaction(enzymes, spont, rid, cofactors, side, attrs, direction, line=line)


def _parse_body(body: str, line: int, col: int) -> list[Reaction]:
    """The reactions of a link; an empty body is one reaction without ID and enzyme."""
    return [_parse_reaction(part, line, col) for part in _split_reactions(body)]


# ---------------------------------------------------------------------------------- parser


class _Parser:
    def __init__(self, text: str):
        self.toks, self.comments = _lex(text)
        self.i = 0

    # --- helpers
    def peek(self, k: int = 0) -> _Tok:
        return self.toks[min(self.i + k, len(self.toks) - 1)]

    def error(self, message: str, tok: _Tok | None = None) -> PWLError:
        tok = tok or self.peek()
        return PWLError(message, tok.line, tok.col)

    def take(self, kind: str | None = None, value=None) -> _Tok:
        t = self.peek()
        if (kind and t.kind != kind) or (value is not None and t.value != value):
            want = repr(value) if value is not None else kind
            got = t.value if t.kind != "EOF" else "end of file"
            raise self.error(f"expected {want}, got {got!r}")
        self.i += 1
        return t

    def is_punct(self, value: str, k: int = 0) -> bool:
        t = self.peek(k)
        return t.kind == "PUNCT" and t.value == value

    def is_id(self, value: str | None = None, k: int = 0) -> bool:
        t = self.peek(k)
        return t.kind == "ID" and (value is None or t.value == value)

    def word(self, choices) -> str:
        t = self.take("ID")
        if t.value not in choices:
            raise self.error(f"expected one of {', '.join(choices)}, got {t.value!r}", t)
        return t.value

    def integer(self) -> int:
        t = self.take("NUM")
        if not t.value.isdigit():
            raise self.error(f"expected an integer, got {t.value!r}", t)
        return int(t.value)

    # --- document
    def document(self) -> Document:
        if not self.is_id("pwl"):
            raise self.error("a PWL file starts with 'pwl 2.0;'")
        self.take()
        version = self.take("NUM").value
        if version.split(".")[0] != VERSION.split(".")[0]:
            raise self.error(f"unsupported PWL version {version} (this parser reads {VERSION})")
        self.take("PUNCT", ";")
        self.take("ID", "subsystem")
        sid = self.take("ID").value
        title = self.take("STR").value
        doc = Document(sid, title, version)
        self.take("PUNCT", "{")
        while not self.is_punct("}"):
            self.statement(doc)
        self.take("PUNCT", "}")
        if self.peek().kind != "EOF":
            raise self.error("content after the end of the subsystem")
        self.apply_directions(doc)
        return doc

    def apply_directions(self, doc: Document) -> None:
        """Apply the ``direction`` statements, in their order, to the reactions."""
        index = doc.reaction_index()
        for setting in doc.directions:
            for rid in setting.reactions:
                r = index.get(rid)
                if r is None:
                    raise PWLError(f"direction of unknown reaction {rid}", setting.line)
                if setting.direction is not None and not r.sided:
                    raise PWLError(f"reaction {rid} has no sides, so it has no direction", setting.line)
                r.direction = setting.direction

    def statement(self, doc: Document) -> None:
        t = self.peek()
        if t.kind == "ID":
            if t.value == "flow":
                self.take()
                if doc.flow is not None:
                    raise self.error("flow given twice", t)
                doc.flow = self.word(FLOWS)
                self.take("PUNCT", ";")
                return
            if t.value == "compartment":
                self.take()
                name = self.take("ID").value
                label = self.take("STR").value if self.peek().kind == "STR" else None
                if name in doc.compartments:
                    raise self.error(f"compartment {name!r} declared twice", t)
                doc.compartments[name] = label
                self.take("PUNCT", ";")
                return
            if t.value == "compound":
                self.take()
                token = self.take("ID").value
                if token in KEYWORDS:
                    raise self.error(f"{token!r} is a keyword and cannot be a compound", t)
                if token in doc.compounds:
                    raise self.error(f"compound {token!r} declared twice", t)
                label = self.take("STR").value if self.peek().kind == "STR" else None
                xrefs = []
                while self.is_id() and self.is_punct(":", 1):
                    ns = self.take().value
                    self.take()
                    v = self.take()
                    if v.kind not in ("ID", "NUM"):
                        raise self.error("expected an identifier after ':'", v)
                    xrefs.append(f"{ns}:{v.value}")
                doc.compounds[token] = Compound(token, label, xrefs, t.line)
                self.take("PUNCT", ";")
                return
            if t.value == "port":
                self.take()
                p = self.participant(allow_coef=False)
                self.take("ID", "side")
                side = self.word(PORT_SIDES)
                target = None
                if self.is_id("to"):
                    self.take()
                    target = self.take("ID").value
                self.take("PUNCT", ";")
                doc.ports.append(Port(p, side, target, t.line))
                return
            if t.value == "ring":
                self.take()
                sense = self.word(SENSES)
                self.take("PUNCT", "{")
                body = []
                while not self.is_punct("}"):
                    body.append(self.ring_statement())
                self.take("PUNCT", "}")
                if not body:
                    raise self.error("empty ring", t)
                doc.statements.append(Ring(sense, body, t.line))
                return
            if t.value == "star":
                doc.statements.append(self.star())
                return
            if t.value == "react":
                doc.statements.append(self.react())
                return
            if t.value == "motif":
                doc.motifs.append(self.motif())
                return
            if t.value == "direction":
                self.take()
                rids = [self.take("ID").value]
                while self.is_punct(","):
                    self.take()
                    rids.append(self.take("ID").value)
                if self.is_punct("?"):
                    self.take()
                    direction = None
                else:
                    direction = self.take("DIR").value
                self.take("PUNCT", ";")
                doc.directions.append(DirectionSetting(rids, direction, t.line))
                return
            if t.value == "effector":
                self.take()
                p = self.participant(allow_coef=False)
                effect = self.word(EFFECTS)
                rids = [self.take("ID").value]
                while self.is_punct(","):
                    self.take()
                    rids.append(self.take("ID").value)
                self.take("PUNCT", ";")
                doc.effectors.append(Effector(p, effect, rids, t.line))
                return
            if t.value in ("chord", "pwl", "subsystem", "wrap"):
                raise self.error(f"{t.value!r} is not allowed here", t)
        doc.statements.append(self.chain_statement())

    def ring_statement(self):
        if self.is_id("chord"):
            self.take()
            chain = self.chain()
            chain.chord = True
            for link in chain.links:
                for r in link.reactions:
                    r.context = "chord"
            self.take("PUNCT", ";")
            return chain
        if self.is_id("react"):
            return self.react()
        if self.is_id() and self.peek().value in KEYWORDS:
            raise self.error(f"{self.peek().value!r} is not allowed inside a ring")
        chain = self.chain("ring")
        self.take("PUNCT", ";")
        return chain

    def chain_statement(self) -> Chain:
        chain = self.chain()
        if self.is_id("wrap"):
            self.take()
            kind = self.word(WRAPS)
            self.take("ID", "period")
            chain.wrap = (kind, self.integer())
        self.take("PUNCT", ";")
        return chain

    def star(self) -> Star:
        t = self.take("ID", "star")
        hub = self.node()
        if hub.open or hub.group:
            raise self.error("the hub of a star is a single compound", t)
        self.take("ID", "arc")
        arc = self.integer()
        self.take("ID", "order")
        order = self.word(SENSES)
        self.take("PUNCT", "{")
        spokes = []
        while not self.is_punct("}"):
            line = self.peek().line
            nodes, links = [hub], []
            while True:
                links.append(self.link("star"))
                nodes.append(self.node())
                if not self.peek().kind == "BRK":
                    break
            self.take("PUNCT", ";")
            spokes.append(Chain(nodes, links, line=line))
            self._connect(spokes[-1])
        self.take("PUNCT", "}")
        return Star(hub, arc, order, spokes, t.line)

    def react(self) -> React:
        t = self.take("ID", "react")
        b = self.take("BRK")
        pre, turn, body, post = b.value
        if pre or turn or post:
            raise self.error("a react statement takes a plain [...] without arrows", b)
        reactions = _parse_body(body, b.line, b.col)
        if len(reactions) != 1:
            raise self.error("a react statement has one reaction", b)
        r = reactions[0]
        if r.direction:
            raise self.error("the participants of a react statement have no sides, so it has no direction", b)
        r.context = "react"
        r.comment = self.comments.get(b.line, "")
        self.take("PUNCT", "{")
        participants = [self.participant()]
        while self.is_punct(","):
            self.take()
            participants.append(self.participant())
        self.take("PUNCT", "}")
        self.take("PUNCT", ";")
        if len(participants) < 2:
            raise self.error("a react statement lists at least two participants", t)
        r.unsided = participants
        return React(r, t.line)

    def motif(self) -> Motif:
        t = self.take("ID", "motif")
        kind = self.take("ID").value
        params = {}
        while self.is_id() and self.is_punct("=", 1):
            key = self.take().value
            self.take()
            v = self.take()
            if v.kind not in ("ID", "NUM", "STR"):
                raise self.error("expected a value", v)
            params[key] = v.value
        members = self.member_list()
        roles = {}
        while self.is_id():
            role = self.take().value
            if role in roles:
                raise self.error(f"role {role!r} given twice", t)
            roles[role] = self.member_list()
        self.take("PUNCT", ";")
        return Motif(kind, members, params, roles, t.line)

    def member_list(self) -> list[str]:
        self.take("PUNCT", "{")
        out = []
        while not self.is_punct("}"):
            tok = self.take("ID").value
            if self.is_punct("@"):
                self.take()
                tok += "@" + self.take("ID").value
            out.append(tok)
        self.take("PUNCT", "}")
        return out

    # --- chains
    def chain(self, context: str = "chain") -> Chain:
        line = self.peek().line
        nodes = [self.node()]
        links = []
        while self.peek().kind == "BRK":
            links.append(self.link(context))
            nodes.append(self.node())
        if not links:
            raise self.error("expected a link -[...]- after the compound")
        chain = Chain(nodes, links, line=line)
        self._connect(chain)
        return chain

    @staticmethod
    def _connect(chain: Chain) -> None:
        for k, link in enumerate(chain.links):
            for r in link.reactions:
                r.left = list(chain.nodes[k].participants)
                r.right = list(chain.nodes[k + 1].participants)

    def link(self, context: str) -> Link:
        t = self.take("BRK")
        pre, turn, body, post = t.value
        if pre is None or post is None:
            raise self.error("a link is written -[...]- (with arrows -[...]->, <-[...]-, <-[...]->, -[...]=>, "
                             "<=[...]-)", t)
        if (pre, post) not in LINK_ARROWS:
            raise self.error(f"undefined link arrows {pre}[...]{post}", t)
        arrow = LINK_ARROWS[(pre, post)]
        reactions = _parse_body(body, t.line, t.col)
        comment = self.comments.get(t.line, "")
        for r in reactions:
            # a direction of the reaction itself overrides the arrows of the link
            if r.direction is None:
                r.direction = arrow
            r.comment = comment
            r.context = context
        return Link(reactions, turn, arrow, t.line)

    def node(self) -> Node:
        if self.is_id("_"):
            self.take()
            return Node([])
        if self.is_punct("("):
            t = self.take()
            parts = [self.participant()]
            while self.is_punct("+"):
                self.take()
                parts.append(self.participant())
            self.take("PUNCT", ")")
            if len(parts) < 2:
                raise self.error("a group has at least two members", t)
            return Node(parts, group=True)
        return Node([self.participant()])

    def participant(self, allow_coef: bool = True) -> Participant:
        coef = None
        t = self.peek()
        if t.kind in ("NUM", "ID") and self.is_punct("*", 1):
            if not allow_coef:
                raise self.error("no coefficient allowed here")
            coef = self.take().value
            self.take()
            if not (_NUM.fullmatch(coef) or coef == "n"):
                raise self.error(f"bad coefficient {coef!r}", t)
        t = self.take("ID")
        if t.value in KEYWORDS:
            raise self.error(f"keyword {t.value!r} where a compound was expected", t)
        if t.value == "_":
            raise self.error("'_' cannot be part of a group or have a coefficient", t)
        copy = 1
        if self.is_punct("#"):
            self.take()
            copy = self.integer()
            if copy < 2:
                raise self.error("copy numbers start at #2", t)
        loc = None
        if self.is_punct("@"):
            self.take()
            loc = self.take("ID").value
        return Participant(t.value, copy, coef, loc)


def parse(text: str, validate: bool = True) -> Document:
    """Parse a PWL 2.0 text; with `validate`, semantic errors raise :class:`PWLError` as well."""
    doc = _Parser(text).document()
    if validate:
        problems = doc.problems()
        if problems:
            more = f" (and {len(problems) - 10} more)" if len(problems) > 10 else ""
            raise PWLError("; ".join(problems[:10]) + more)
    return doc


def parse_file(path, validate: bool = True) -> Document:
    with open(path, encoding="utf-8") as f:
        return parse(f.read(), validate)


# ---------------------------------------------------------------------------------- validation


def validate(doc: Document) -> list[str]:
    """Semantic errors of a document.

    - reaction IDs are unique (reactions need no ID, compounds and compartments need no declaration)
    - motifs and effectors refer to existing reactions/compounds; members of the known motif kinds have the right type
    - ports name drawn nodes
    """
    problems = []
    reactions = doc.reactions()
    rids: dict[str, Reaction] = {}
    for r in reactions:
        where = f"line {r.line}"
        if r.rid and r.rid in rids:
            problems.append(f"{where}: reaction ID {r.rid} used twice")
        elif r.rid:
            rids[r.rid] = r
        if r.context != "react" and not r.left and not r.right:
            problems.append(f"{where}: reaction {r.rid or r.enzyme_label} has no main metabolite")
        for value in r.attrs.values():
            if not value:
                problems.append(f"{where}: empty attribute")
    drawn = doc.nodes()
    for port in doc.ports:
        if port.participant.node not in drawn:
            problems.append(f"line {port.line}: port {port.participant.node!r} is not a drawn node")
    species = {p.species for r in reactions for p in r.participants + r.consumed + r.produced}
    tokens = set(doc.compounds)
    for motif in doc.motifs:
        where = f"line {motif.line}"
        spec = MOTIF_KINDS.get(motif.kind)
        if spec is None:
            continue  # free motif kinds are allowed, their members are not checked
        member_type, role_types = spec
        if not motif.members:
            problems.append(f"{where}: motif {motif.kind} without members")
        for role in motif.roles:
            if role not in role_types:
                problems.append(f"{where}: motif {motif.kind} has no role {role!r}")
        for role, members in [("", motif.members)] + list(motif.roles.items()):
            kind = member_type if not role else role_types.get(role)
            for m in members:
                if kind == "reaction" and m not in rids:
                    problems.append(f"{where}: motif {motif.kind}: unknown reaction {m}")
                if kind == "compound" and m not in species and m.split("@")[0] not in tokens:
                    problems.append(f"{where}: motif {motif.kind}: unknown compound {m}")
    for eff in doc.effectors:
        for rid in eff.reactions:
            if rid not in rids:
                problems.append(f"line {eff.line}: effector of unknown reaction {rid}")
    return problems
