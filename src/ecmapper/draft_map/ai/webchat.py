"""Web chats as AI connections: the tool loop of an AI job through the text of a chat on a web page.

A web chat has no tool calls of its own: ecmapper's first message tells the
model the tools and how to call them (:data:`PROTOCOL`), the model writes its
calls as JSON blocks in its reply (:func:`parse`) and gets their results in
ecmapper's next message. The messages are numbered: message N ends with the
line :func:`tail` and the model ends its reply with the line :func:`end_line`,
so a transport that reads the page knows where the reply is and when it is
complete (:func:`reply_of`). A transport carries the messages:
``exchange(message, number, images, stop)`` returns the reply, with its end
line if it is complete (None: stopped),
and has `images` (True: it can send images);
the GUI has one that drives a browser and one for copying by hand
(:mod:`..gui.webchat`).

The browser may not be used for chats whose terms forbid automated use
(:data:`BLOCKED`); copying by hand is fine there.
"""

from __future__ import annotations

import json
import re
from urllib.parse import urlparse

from .providers import Chat, Connection, ProviderError, ToolCall, ToolResult, Turn

#: the transports of a web chat connection (:attr:`.providers.Connection.transport`)
BROWSER, CLIPBOARD = "browser", "clipboard"

_TERMS = ("Use an API connection or a coding agent's work folder (AI menu) with your subscription instead, or the "
          "web chat with copy/paste.")
#: domain: why the browser may not drive chats there
BLOCKED = {
    "claude.ai": "Anthropic's terms do not allow automated use of claude.ai (only with an API key)",
    "claude.com": "Anthropic's terms do not allow automated use of Claude's apps (only with an API key)",
    "chatgpt.com": "OpenAI's terms do not allow extracting the output of ChatGPT programmatically",
    "chat.openai.com": "OpenAI's terms do not allow extracting the output of ChatGPT programmatically",
    "gemini.google.com": "Google's terms do not allow automated use of Gemini",
    "aistudio.google.com": "Google's terms do not allow automated use of AI Studio's chat",
}

PROTOCOL = """You are in a conversation with ecmapper, a program on the user's computer: ecmapper writes the \
messages (each starts with "ecmapper message N") and reads your replies; the user's requests come through it. You \
need no files and run no code yourself: the map and the toolkit are on the user's computer, and ecmapper runs the \
tools below there.

{system}

=== How to call the tools ===
ecmapper runs these tools for you and sends their results in its next message:

{tools}

Write each call as a JSON object in a fenced code block, for example

```json
{{"tool": "list_files", "input": {{"path": "."}}}}
```

You may write several calls in one reply: they run in this order and all results come in the next message. Use only \
these tools: do not run code yourself (no code interpreter, analysis tool or canvas) and never make up results. Write \
no other JSON objects with a "tool" key. A reply without calls ends your turn, and the user reads it.

Chats cut long replies off: keep each reply short (at most about 150 lines of a file per reply) and write a longer \
file in parts, the first with write_file, each further one with write_file and "append": true in a later reply \
(only write_file has "append"; edit_file changes a part of a file).

End every reply with the line @@END N@@, N being the number of the message you answer."""


def tail(number: int) -> str:
    """The last line of ecmapper's message `number`."""
    return f"(end of ecmapper message {number}: reply below)"


def end_line(number: int) -> str:
    """The line that ends the reply to message `number`."""
    return f"@@END {number}@@"


def blocked(url: str) -> str | None:
    """Why the browser may not drive the chat at `url` (None: it may)."""
    host = (urlparse(url if "//" in url else f"https://{url}").hostname or "").lower()
    for domain, why in BLOCKED.items():
        if host == domain or host.endswith("." + domain):
            return f"{why}. {_TERMS}"
    return None


def reply_of(page: str, number: int) -> tuple[str, bool]:
    """The reply to message `number` in the text of a chat (the page, or what the user copied) and whether it is
    complete (ends with :func:`end_line`). It starts after the message's :func:`tail` (or after the end of the reply
    before, if the chat shows the message as an attachment)."""
    start = page.rfind(tail(number))
    if start >= 0:
        start += len(tail(number))
    else:
        start = page.rfind(end_line(number - 1))
        start = start + len(end_line(number - 1)) if start >= 0 else 0
    end = page.find(end_line(number), start)
    return (page[start:end] if end >= 0 else page[start:]).strip(), end >= 0


_FENCE = re.compile(r"^[ \t]*```[\w-]*[ \t]*$", re.MULTILINE)
_NAME = re.compile(r'"tool"\s*:\s*"([^"]*)"')
_CALL = re.compile(r'\{\s*"tool"\s*:')


CUT_OFF = ("your reply was cut off in this call (the chat's limit of the length of a reply), so it was not run: keep "
           "replies short and write a long file in parts (write_file, then write_file with \"append\": true, one part "
           "per reply)")


def parse(reply: str, number: int, complete: bool = True) -> tuple[str, list[ToolCall]]:
    """The text of a reply (without its calls) and the calls in it (IDs N.1, N.2, …). A call that is not valid JSON
    or has no object as input is a call with an error (the agent tells the model); the last one of a reply that is
    not `complete` was cut off (:data:`CUT_OFF`)."""
    decoder = json.JSONDecoder(strict=False)  # line breaks and tabs in strings (models write files so)
    calls, text, i, last = [], [], 0, 0
    while True:
        i = reply.find("{", i)
        if i < 0:
            break
        try:
            value, end = decoder.raw_decode(reply, i)
        except json.JSONDecodeError as e:
            value, end, problem = None, -1, e
            p = e.pos
            if complete and e.msg.startswith("Expecting ',' delimiter") and re.match(r"\s*($|\{|`)", reply[p:]):
                # a call without its last brace(s), before the next call (models miscount them when the content is JSON)
                for k in (1, 2):
                    try:
                        value, end = decoder.raw_decode(reply[:p] + "}" * k + reply[p:], i)
                        end -= k
                        break
                    except json.JSONDecodeError:
                        pass
        if isinstance(value, dict) and "tool" in value:
            name, data = value.get("tool"), value.get("input", {})
            error = None if isinstance(data, dict) else "'input' must be a JSON object"
            error = error if isinstance(name, str) else "'tool' must be the name of a tool"
            calls.append(ToolCall(f"{number}.{len(calls) + 1}", name if isinstance(name, str) else "?",
                                  data if isinstance(data, dict) else {}, error))
        elif value is None and _CALL.match(reply, i):
            match = _NAME.search(reply, i, i + 200)
            line, near = reply.count("\n", i, problem.pos) + 1, reply[problem.pos:problem.pos + 30]
            where = f"line {line} of the call, at {near!r}"
            calls.append(ToolCall(f"{number}.{len(calls) + 1}", match.group(1) if match else "?", {},
                                  f"the call is not valid JSON ({problem.msg} in {where}): write it as one JSON "
                                  "object"))
            following = _CALL.search(reply, problem.pos)  # the next call (chats show code without its fences)
            end = reply.find("\n\n", i)
            end = min(len(reply) if end < 0 else end, following.start() if following else len(reply))
        else:
            i += 1
            continue
        text.append(reply[last:i])
        i = last = end
    text.append(reply[last:])
    if not complete and calls and calls[-1].error:
        calls[-1].error = CUT_OFF
    shown = _FENCE.sub("", "".join(text))
    return re.sub(r"\n{3,}", "\n\n", shown).strip(), calls


class WebChat(Chat):
    """A conversation in a web chat (see the module); `transport` carries the messages. `history`: the messages and
    replies (the chat itself keeps the conversation)."""

    def __init__(self, connection: Connection, system: str, tools, transport):
        super().__init__(connection, system, tools)
        self.transport = transport
        self.number = 0
        self.pending: list[str] = []
        self.images: list[bytes] = []
        self.names: dict[str, str] = {}  # call ID -> tool
        self.history: list[tuple[str, str]] = []

    def intro(self) -> str:
        tools = "\n".join(f"- {t.name}: {t.description}\n  input (JSON schema): {json.dumps(t.schema)}"
                          for t in self.tool_list)
        return PROTOCOL.format(system=self.system, tools=tools)

    def add_user(self, text: str) -> None:
        self.pending.append(f"=== The user ===\n{text}")

    def add_results(self, results: list[ToolResult]) -> None:
        parts, left_out = [], 0
        images = self.connection.images and getattr(self.transport, "images", False)
        for r in results:
            head = f"[call {r.id}: {self.names.get(r.id, 'tool')}{' - ERROR' if r.is_error else ''}]"
            parts.append(f"{head}\n{r.text or '(no output)'}")
            if images:
                self.images += r.images
            else:
                left_out += len(r.images)
        note = f"\n\n({left_out} image(s) of the results can not be shown in this chat.)" if left_out else ""
        self.pending.append("=== Results of your calls ===\n" + "\n\n".join(parts) + note)

    def step(self, on_text=None, on_thinking=None, stop=None) -> Turn:
        self.number += 1
        n = self.number
        body = "\n\n".join(self.pending)
        message = f"ecmapper message {n}\n\n" + (f"{self.intro()}\n\n" if n == 1 else "") + f"{body}\n\n{tail(n)}"
        images, self.pending, self.images = self.images, [], []
        reply = self.transport.exchange(message, n, images, stop)
        if reply is None:
            return Turn("", [], "cancelled")
        complete = end_line(n) in reply
        reply = reply.replace(end_line(n), "").strip()
        self.history.append((message, reply))
        text, calls = parse(reply, n, complete)
        self.names = {c.id: c.name for c in calls}
        if text and on_text is not None:
            on_text(text)
        broken = any(c.error for c in calls)  # the reply as it came, for the transcript
        return Turn(text, calls, "tool_use" if calls else "end", detail=reply if broken else "")


def check_connection(connection: Connection) -> None:
    """ProviderError if the web chat `connection` can not be used (a blocked chat in the browser)."""
    if connection.transport == BROWSER:
        if not connection.base_url:
            raise ProviderError("Give the URL of the web chat (its page) for the browser.")
        reason = blocked(connection.base_url)
        if reason:
            raise ProviderError(reason)
