"""Language-model back-ends. All of them just turn a list of chat messages into reply text.

  OpenAICompat  any OpenAI-style /chat/completions endpoint (OpenAI, OpenRouter, Azure-compatible gateways, Ollama,
                llama.cpp, vLLM, LM Studio, ...). Only the standard library is used.
  Anthropic     the Anthropic Messages API.
  Manual        no API at all: writes each prompt to a file, you paste it into ANY chat UI and paste the reply back.
  Scripted      replays prepared replies (tests / demos).
The harness never gives a model tools, file access or code execution; replies are parsed as JSON only.
"""
import copy, json, os, sys, time, urllib.request, urllib.error
from urllib.parse import urlparse
from .._io import read_text, write_text


class LLMError(Exception):
    pass


class LLM:
    name = 'llm'
    def chat(self, messages, turn=0):
        raise NotImplementedError


def _post(url, headers, body, timeout, retries=3):
    if urlparse(url).scheme not in ('http', 'https'): raise LLMError('only http(s) endpoints are allowed')
    data = json.dumps(body).encode()
    last = None
    for attempt in range(retries):
        req = urllib.request.Request(url, data=data, headers={**headers, 'Content-Type': 'application/json'}, method='POST')
        try:
            with urllib.request.urlopen(req, timeout=timeout) as r:
                return json.loads(r.read().decode('utf-8'))
        except urllib.error.HTTPError as e:
            last = f'HTTP {e.code}: {e.read().decode("utf-8", "ignore")[:300]}'
            if e.code not in (408, 429, 500, 502, 503, 504): break
        except (urllib.error.URLError, TimeoutError, ConnectionError) as e:
            last = f'connection problem: {e}'
        time.sleep(2 * (attempt + 1))
    raise LLMError(f'request failed ({last})')


class OpenAICompat(LLM):
    def __init__(self, base_url, model, api_key=None, temperature=0.2, max_tokens=8000, json_mode=False, timeout=300):
        self.url = base_url.rstrip('/') + '/chat/completions'; self.model = model; self.key = api_key
        self.temperature = temperature; self.max_tokens = max_tokens; self.json_mode = json_mode; self.timeout = timeout
        self.name = f'openai-compatible:{model}'
    def chat(self, messages, turn=0):
        body = {'model': self.model, 'messages': messages, 'temperature': self.temperature, 'max_tokens': self.max_tokens}
        if self.json_mode: body['response_format'] = {'type': 'json_object'}
        headers = {'Authorization': f'Bearer {self.key}'} if self.key else {}
        out = _post(self.url, headers, body, self.timeout)
        try: return out['choices'][0]['message']['content'] or ''
        except (KeyError, IndexError, TypeError): raise LLMError('unexpected response format from the endpoint')


class Anthropic(LLM):
    def __init__(self, model, api_key, base_url='https://api.anthropic.com', temperature=0.2, max_tokens=8000, timeout=300):
        self.url = base_url.rstrip('/') + '/v1/messages'; self.model = model; self.key = api_key
        self.temperature = temperature; self.max_tokens = max_tokens; self.timeout = timeout; self.name = f'anthropic:{model}'
    def chat(self, messages, turn=0):
        system = '\n\n'.join(m['content'] for m in messages if m['role'] == 'system')
        msgs = [m for m in messages if m['role'] != 'system']
        body = {'model': self.model, 'max_tokens': self.max_tokens, 'temperature': self.temperature, 'messages': msgs}
        if system: body['system'] = system
        out = _post(self.url, {'x-api-key': self.key or '', 'anthropic-version': '2023-06-01'}, body, self.timeout)
        try: return ''.join(b.get('text', '') for b in out['content'] if b.get('type') == 'text')
        except (KeyError, TypeError): raise LLMError('unexpected response format from the endpoint')


class Scripted(LLM):
    """Replays prepared replies in order (str or dict). Raises when the script is exhausted."""
    name = 'scripted'
    def __init__(self, replies): self.replies = list(replies); self.i = 0; self.seen = []
    def chat(self, messages, turn=0):
        self.seen.append(copy.deepcopy(messages))
        if self.i >= len(self.replies): raise LLMError('scripted LLM ran out of replies')
        r = self.replies[self.i]; self.i += 1
        return r if isinstance(r, str) else json.dumps(r)


class Manual(LLM):
    """Copy/paste mode: each prompt goes to <dir>/turn_NNN_prompt.md; the reply is read from turn_NNN_reply.txt
    (if it exists) or pasted into the terminal and finished with a line containing only END."""
    name = 'manual'
    def __init__(self, directory, stdin=None, out=None):
        self.dir = directory; os.makedirs(directory, exist_ok=True); self.stdin = stdin or sys.stdin; self.out = out or sys.stdout
    def chat(self, messages, turn=0):
        pf = os.path.join(self.dir, f'turn_{turn:03d}_prompt.md'); rf = os.path.join(self.dir, f'turn_{turn:03d}_reply.txt')
        parts = []
        for m in messages:
            head = {'system': '# SYSTEM INSTRUCTIONS', 'user': '# USER MESSAGE', 'assistant': '# YOUR PREVIOUS REPLY'}[m['role']]
            parts.append(f"{head}\n\n{m['content']}")
        write_text(pf, '\n\n---\n\n'.join(parts) + '\n')
        if os.path.exists(rf): return read_text(rf)
        print(f'\n=== Turn {turn}: open {pf}, paste its whole content into your chat model, then paste the model\'s reply here.\n'
              f'    (finish with a line containing only END; or save the reply as {rf} and press Enter)', file=self.out, flush=True)
        lines = []
        for line in self.stdin:
            if line.strip() == 'END': break
            lines.append(line)
            if len(lines) == 1 and line.strip() == '' and os.path.exists(rf): return read_text(rf)
        text = ''.join(lines)
        if not text.strip() and os.path.exists(rf): return read_text(rf)
        write_text(rf, text)
        return text
