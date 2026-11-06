"""The tool loop of an AI job: the model answers, the tools it calls run (after approval where needed), their
results go back in one message, until the model calls no more tools, the user stops it or `max_steps` is reached."""

from __future__ import annotations

import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from .providers import Chat, ToolCall, ToolResult, Turn


class ToolError(Exception):
    """A problem the model can fix: the message is the (error) result of the tool."""


@dataclass
class Tool:
    """A tool for the model: `run(**input)` returns a text or (text, [PNG images]); `needs_approval`: the user
    approves each call (see :meth:`Events.approve`)."""

    name: str
    description: str
    schema: dict
    run: Callable
    needs_approval: bool = False


class Events:
    """What the agent reports and asks while it runs (the GUI overrides these; the defaults do nothing)."""

    def text(self, delta: str) -> None:
        pass

    def thinking(self, delta: str) -> None:
        pass

    def turn(self, turn: Turn) -> None:
        pass

    def tool_started(self, call: ToolCall) -> None:
        pass

    def tool_finished(self, call: ToolCall, result: ToolResult) -> None:
        pass

    def approve(self, call: ToolCall) -> bool:
        return True

    def cancelled(self) -> bool:
        return False


TYPES = {"string": str, "integer": int, "number": (int, float), "boolean": bool, "array": list, "object": dict}


def schema_problem(value, schema: dict, where: str = "input") -> str | None:
    """Why `value` does not fit the JSON schema (types, required and unknown properties, items), or None."""
    kind = schema.get("type")
    if kind in TYPES and (not isinstance(value, TYPES[kind]) or (kind in ("integer", "number")
                                                                  and isinstance(value, bool))):
        return f"{where} must be of type {kind}"
    if kind == "object":
        props = schema.get("properties", {})
        for name in schema.get("required", []):
            if name not in value:
                return f"{where}: {name!r} is missing"
        if schema.get("additionalProperties") is False:
            extra = sorted(set(value) - set(props))
            if extra:
                return f"{where}: unknown {', '.join(map(repr, extra))}"
        for name, sub in props.items():
            if name in value:
                problem = schema_problem(value[name], sub, f"{where}.{name}")
                if problem:
                    return problem
    if kind == "array" and "items" in schema:
        for i, item in enumerate(value):
            problem = schema_problem(item, schema["items"], f"{where}[{i}]")
            if problem:
                return problem
    return None


class Agent:
    """Runs the tools of a chat; `transcript`: a JSON lines file that records the conversation."""

    def __init__(self, chat: Chat, tools: list[Tool], transcript: Path | None = None, max_steps: int = 80):
        self.chat, self.max_steps = chat, max_steps
        self.tools = {t.name: t for t in tools}
        self.transcript = transcript

    def log(self, **entry) -> None:
        if self.transcript is not None:
            with open(self.transcript, "a", encoding="utf-8") as f:
                f.write(json.dumps(dict(time=round(time.time(), 1), **entry), default=str) + "\n")

    def run(self, text: str | None, events: Events | None = None) -> str:
        """Add the user's `text` (None: continue) and run until the model is done. Returns why it stopped: done,
        stopped (by the user), max_steps, max_tokens or refusal (with the category: "refusal: cyber")."""
        events = events or Events()
        if text is not None:
            self.chat.add_user(text)
            self.log(role="user", text=text)
        for _ in range(self.max_steps):
            if events.cancelled():
                return "stopped"
            turn = self.chat.step(on_text=events.text, on_thinking=events.thinking, stop=events.cancelled)
            if turn.stop == "cancelled":
                return "stopped"
            self.log(role="assistant", text=turn.text, thinking=turn.thinking, stop=turn.stop, usage=turn.usage,
                     calls=[dict(id=c.id, name=c.name, input=c.input, error=c.error) for c in turn.calls])
            events.turn(turn)
            if turn.stop == "refusal":
                return f"refusal: {turn.detail}" if turn.detail else "refusal"
            if turn.stop == "pause":
                continue
            if not turn.calls:
                return "max_tokens" if turn.stop == "max_tokens" else "done"
            results = []
            for call in turn.calls:
                if events.cancelled():
                    results.append(ToolResult(call.id, "Stopped by the user before this call ran.", is_error=True))
                    continue
                results.append(self.execute(call, events))
            self.chat.add_results(results)
            self.log(role="tool", results=[dict(id=r.id, text=r.text, images=len(r.images), error=r.is_error)
                                           for r in results])
            if events.cancelled():
                return "stopped"
        return "max_steps"

    def execute(self, call: ToolCall, events: Events) -> ToolResult:
        tool = self.tools.get(call.name)
        if tool is None:
            return ToolResult(call.id, f"Unknown tool {call.name!r}; the tools are {', '.join(self.tools)}.",
                              is_error=True)
        problem = call.error or schema_problem(call.input, tool.schema)
        if problem:
            return ToolResult(call.id, f"Invalid call of {call.name}: {problem}. Call it again.", is_error=True)
        if tool.needs_approval and not events.approve(call):
            result = ToolResult(call.id, "The user declined this call.", is_error=True)
            events.tool_finished(call, result)
            return result
        events.tool_started(call)
        try:
            out = tool.run(**call.input)
            text, images = (out if isinstance(out, tuple) else (out, []))
            result = ToolResult(call.id, str(text), list(images))
        except ToolError as e:
            result = ToolResult(call.id, str(e), is_error=True)
        except Exception as e:  # noqa: BLE001 (a failing tool is reported to the model, not raised)
            result = ToolResult(call.id, f"{type(e).__name__}: {e}", is_error=True)
        events.tool_finished(call, result)
        return result
