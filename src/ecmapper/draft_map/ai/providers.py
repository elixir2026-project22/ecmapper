"""Connections to AI models with tool use: the Anthropic API (anthropic SDK) and OpenAI-compatible APIs (openai SDK,
chat completions: OpenAI, OpenRouter, Mistral, Gemini, local servers such as Ollama and LM Studio); web chats are in
:mod:`.webchat`.

A chat keeps its history in the provider's own format (append-only, as sent) and takes one step at a time
(:meth:`Chat.step`): the model's answer with the tools it calls; the caller runs them and adds the results
(:meth:`Chat.add_results`), see :mod:`.agent`.
"""

from __future__ import annotations

import base64
import json
import os
from dataclasses import dataclass, field

#: the kinds of connections
ANTHROPIC, OPENAI, WEBCHAT = "anthropic", "openai", "webchat"
#: Anthropic models offered in the GUI (the first is the default)
ANTHROPIC_MODELS = ["claude-opus-5-5", "claude-sonnet-5-5", "claude-haiku-5-5"]
#: Anthropic models with server-side fallbacks after a refusal (fallbacks="default")
FALLBACK_MODELS = {"claude-fable-5-1", "claude-opus-5-5", "claude-opus-5", "claude-sonnet-5-5"}
FALLBACK_BETA = "server-side-fallback-2026-07-01"
EFFORTS = ["low", "medium", "high", "xhigh", "max"]
#: presets of connections: name -> (kind, base URL, environment variable of the key, model)
PRESETS = {
    "Anthropic": (ANTHROPIC, "", "ANTHROPIC_API_KEY", ANTHROPIC_MODELS[0]),
    "OpenAI": (OPENAI, "https://api.openai.com/v1", "OPENAI_API_KEY", ""),
    "OpenRouter": (OPENAI, "https://openrouter.ai/api/v1", "OPENROUTER_API_KEY", ""),
    "Mistral": (OPENAI, "https://api.mistral.ai/v1", "MISTRAL_API_KEY", ""),
    "Gemini (OpenAI-compatible)": (OPENAI, "https://generativelanguage.googleapis.com/v1beta/openai/",
                                   "GEMINI_API_KEY", ""),
    "Ollama (local)": (OPENAI, "http://localhost:11434/v1", "", ""),
    "LM Studio (local)": (OPENAI, "http://localhost:1234/v1", "", ""),
    "Other OpenAI-compatible": (OPENAI, "", "", ""),
}


@dataclass
class Connection:
    """How to reach a model. `key`: the API key itself (else from the environment variable `key_env`; Anthropic:
    else the SDK's own credentials); `images`: send images (renders) to the model. A web chat (kind
    :data:`WEBCHAT`) has its page as `base_url` and `transport` (browser or clipboard, see :mod:`.webchat`)."""

    name: str
    kind: str = ANTHROPIC
    model: str = ANTHROPIC_MODELS[0]
    base_url: str = ""
    key_env: str = "ANTHROPIC_API_KEY"
    key: str = ""
    effort: str = "high"
    images: bool = True
    max_tokens: int = 64000
    transport: str = "browser"

    def api_key(self) -> str | None:
        return self.key or (os.environ.get(self.key_env) if self.key_env else None) or None

    def settings(self) -> dict:
        """What is saved (no key)."""
        return dict(name=self.name, kind=self.kind, model=self.model, base_url=self.base_url, key_env=self.key_env,
                    effort=self.effort, images=self.images, max_tokens=self.max_tokens, transport=self.transport)


@dataclass
class ToolCall:
    id: str
    name: str
    input: dict
    error: str | None = None  # the input could not be read (e.g. cut off at the output limit)


@dataclass
class ToolResult:
    id: str
    text: str
    images: list[bytes] = field(default_factory=list)  # PNG
    is_error: bool = False


@dataclass
class Turn:
    """An answer of the model. `stop`: tool_use, end, max_tokens, refusal, pause (continue with another step) or
    cancelled (not added to the history); `detail`: e.g. the category of a refusal."""

    text: str
    calls: list[ToolCall]
    stop: str
    usage: dict = field(default_factory=dict)
    thinking: str = ""
    detail: str = ""


class ProviderError(Exception):
    """The connection failed (missing SDK, credentials, network, API error); the message says why."""


class Chat:
    """A conversation with a model that can call `tools` (objects with name, description and schema)."""

    def __init__(self, connection: Connection, system: str, tools):
        self.connection, self.system, self.tool_list = connection, system, list(tools)

    def add_user(self, text: str) -> None:
        raise NotImplementedError

    def add_results(self, results: list[ToolResult]) -> None:
        raise NotImplementedError

    def step(self, on_text=None, on_thinking=None, stop=None) -> Turn:
        raise NotImplementedError


def _png(data: bytes) -> str:
    return base64.standard_b64encode(data).decode("ascii")


class AnthropicChat(Chat):
    """Claude through the Messages API (streamed): adaptive thinking (summaries shown), `effort`, prompt caching,
    server-side fallbacks after refusals (models of :data:`FALLBACK_MODELS`); the history is appended unchanged."""

    def __init__(self, connection: Connection, system: str, tools, client=None):
        super().__init__(connection, system, tools)
        if client is None:
            try:
                import anthropic
            except ImportError as e:
                raise ProviderError("the Anthropic connection needs the anthropic package (pip install "
                                    "anthropic)") from e
            client = anthropic.Anthropic(api_key=connection.api_key(), base_url=connection.base_url or None)
        self.client = client
        self.messages: list[dict] = []
        self.tools = [dict(name=t.name, description=t.description, input_schema=t.schema, eager_input_streaming=True)
                      for t in self.tool_list]

    def add_user(self, text: str) -> None:
        self.messages.append({"role": "user", "content": text})

    def add_results(self, results: list[ToolResult]) -> None:
        content = []
        for r in results:
            blocks = [{"type": "text", "text": r.text or "(no output)"}]
            if self.connection.images:
                blocks += [{"type": "image", "source": {"type": "base64", "media_type": "image/png", "data": _png(p)}}
                           for p in r.images]
            content.append({"type": "tool_result", "tool_use_id": r.id, "content": blocks, "is_error": r.is_error})
        self.messages.append({"role": "user", "content": content})

    def request(self) -> dict:
        c = self.connection
        kwargs = dict(model=c.model, max_tokens=c.max_tokens, system=self.system, tools=self.tools,
                      messages=self.messages, thinking={"type": "adaptive", "display": "summarized"},
                      output_config={"effort": c.effort}, cache_control={"type": "ephemeral"})
        if c.model in FALLBACK_MODELS:
            kwargs.update(betas=[FALLBACK_BETA], fallbacks="default")
        return kwargs

    def step(self, on_text=None, on_thinking=None, stop=None) -> Turn:
        try:
            with self.client.beta.messages.stream(**self.request()) as stream:
                for event in stream:
                    if stop is not None and stop():
                        return Turn("", [], "cancelled")
                    if event.type == "content_block_delta":
                        if event.delta.type == "text_delta" and on_text is not None:
                            on_text(event.delta.text)
                        elif event.delta.type == "thinking_delta" and on_thinking is not None:
                            on_thinking(event.delta.thinking)
                message = stream.get_final_message()
        except Exception as e:  # noqa: BLE001 (the SDK's errors: shown to the user)
            raise ProviderError(_api_error(e)) from e
        stop_reason = message.stop_reason
        text = "".join(b.text for b in message.content if b.type == "text")
        thinking = "".join(getattr(b, "thinking", "") or "" for b in message.content if b.type == "thinking")
        usage = {k: getattr(message.usage, k, 0) or 0 for k in
                 ("input_tokens", "output_tokens", "cache_read_input_tokens", "cache_creation_input_tokens")}
        if stop_reason == "refusal":
            details = getattr(message, "stop_details", None)
            return Turn(text, [], "refusal", usage, thinking, getattr(details, "category", None) or "")
        self.messages.append({"role": "assistant", "content": message.content})
        calls = []
        for b in message.content:
            if b.type == "tool_use":
                good = isinstance(b.input, dict)
                error = None if good and stop_reason != "max_tokens" else "the call was cut off (output limit)"
                calls.append(ToolCall(b.id, b.name, b.input if good else {}, error))
        kind = {"tool_use": "tool_use", "end_turn": "end", "max_tokens": "max_tokens", "pause_turn": "pause"}
        return Turn(text, calls, kind.get(stop_reason, "end"), usage, thinking)


class OpenAIChat(Chat):
    """A model of an OpenAI-compatible API (chat completions, streamed, function tools). Images of tool results
    follow as a user message (tool messages take text only)."""

    def __init__(self, connection: Connection, system: str, tools, client=None):
        super().__init__(connection, system, tools)
        if client is None:
            try:
                import openai
            except ImportError as e:
                raise ProviderError("OpenAI-compatible connections need the openai package (pip install openai)") from e
            client = openai.OpenAI(api_key=connection.api_key() or "none", base_url=connection.base_url or None)
        self.client = client
        self.messages: list[dict] = [{"role": "system", "content": system}]
        self.tools = [{"type": "function", "function": dict(name=t.name, description=t.description,
                                                              parameters=t.schema)} for t in self.tool_list]
        self.usage_option = True  # stream_options include_usage (not every server takes it)

    def add_user(self, text: str) -> None:
        self.messages.append({"role": "user", "content": text})

    def add_results(self, results: list[ToolResult]) -> None:
        images = []
        for r in results:
            self.messages.append({"role": "tool", "tool_call_id": r.id,
                                  "content": ("Error: " if r.is_error else "") + (r.text or "(no output)")})
            images += r.images
        if images and self.connection.images:
            content = [{"type": "text", "text": "The images of the tool results above."}]
            content += [{"type": "image_url", "image_url": {"url": f"data:image/png;base64,{_png(p)}"}} for p in images]
            self.messages.append({"role": "user", "content": content})

    def _stream(self):
        kwargs = dict(model=self.connection.model, messages=self.messages, tools=self.tools, stream=True)
        if self.usage_option:
            try:
                return self.client.chat.completions.create(stream_options={"include_usage": True}, **kwargs)
            except Exception as e:  # noqa: BLE001 (a server without stream_options: once more without)
                if "stream_options" not in str(e):
                    raise
                self.usage_option = False
        return self.client.chat.completions.create(**kwargs)

    def step(self, on_text=None, on_thinking=None, stop=None) -> Turn:
        text, thinking, calls, finish, usage = [], [], {}, None, {}
        try:
            for chunk in self._stream():
                if stop is not None and stop():
                    return Turn("", [], "cancelled")
                if getattr(chunk, "usage", None):
                    usage = dict(input_tokens=chunk.usage.prompt_tokens or 0,
                                 output_tokens=chunk.usage.completion_tokens or 0)
                if not chunk.choices:
                    continue
                choice = chunk.choices[0]
                delta = choice.delta
                if delta.content:
                    text.append(delta.content)
                    if on_text is not None:
                        on_text(delta.content)
                reasoning = getattr(delta, "reasoning_content", None) or getattr(delta, "reasoning", None)
                if isinstance(reasoning, str) and reasoning:
                    thinking.append(reasoning)
                    if on_thinking is not None:
                        on_thinking(reasoning)
                for tc in delta.tool_calls or []:
                    call = calls.setdefault(tc.index, dict(id="", name="", arguments=""))
                    call["id"] = tc.id or call["id"]
                    if tc.function is not None:
                        call["name"] += tc.function.name or ""
                        call["arguments"] += tc.function.arguments or ""
                finish = choice.finish_reason or finish
        except Exception as e:  # noqa: BLE001 (the SDK's errors: shown to the user)
            raise ProviderError(_api_error(e)) from e
        message = {"role": "assistant", "content": "".join(text) or None}
        if calls:
            message["tool_calls"] = [{"id": c["id"] or f"call_{i}", "type": "function",
                                      "function": {"name": c["name"], "arguments": c["arguments"] or "{}"}}
                                     for i, c in sorted(calls.items())]
        self.messages.append(message)
        result = []
        for i, c in sorted(calls.items()):
            try:
                args = json.loads(c["arguments"] or "{}")
                error = None if isinstance(args, dict) else "the arguments are not a JSON object"
            except json.JSONDecodeError:
                args, error = {}, "the arguments are not valid JSON" + (" (output limit)" if finish == "length" else "")
            result.append(ToolCall(c["id"] or f"call_{i}", c["name"], args if isinstance(args, dict) else {}, error))
        kind = {"tool_calls": "tool_use", "stop": "end", "length": "max_tokens", "content_filter": "refusal"}
        stop_kind = "tool_use" if result else kind.get(finish, "end")
        return Turn("".join(text), result, stop_kind, usage, "".join(thinking))


def _api_error(e: Exception) -> str:
    """A short message for an error of an SDK."""
    status = getattr(e, "status_code", None)
    text = getattr(e, "message", None) or str(e)
    hints = {401: "the API key is missing or wrong", 403: "the API key lacks permission", 404: "unknown model or URL",
             429: "rate limited, try again later"}
    prefix = f"{type(e).__name__}" + (f" {status}" if status else "")
    hint = hints.get(status, "")
    return f"{prefix}: {text}" + (f" ({hint})" if hint else "")


def make_chat(connection: Connection, system: str, tools, client=None) -> Chat:
    """The chat of the kind of `connection`."""
    if connection.kind == ANTHROPIC:
        return AnthropicChat(connection, system, tools, client)
    if connection.kind == OPENAI:
        return OpenAIChat(connection, system, tools, client)
    if connection.kind == WEBCHAT:
        raise ProviderError("a web chat needs a transport of the GUI (see .webchat.WebChat)")
    raise ProviderError(f"unknown kind of connection: {connection.kind!r}")


def list_models(connection: Connection) -> list[str]:
    """The models the API offers (OpenAI-compatible: /models; Anthropic: the models API)."""
    try:
        if connection.kind == ANTHROPIC:
            import anthropic

            client = anthropic.Anthropic(api_key=connection.api_key(), base_url=connection.base_url or None)
        else:
            import openai

            client = openai.OpenAI(api_key=connection.api_key() or "none", base_url=connection.base_url or None)
        return sorted(m.id for m in client.models.list())
    except ImportError as e:
        raise ProviderError(f"{e.name} is not installed (pip install {e.name})") from e
    except Exception as e:  # noqa: BLE001
        raise ProviderError(_api_error(e)) from e
