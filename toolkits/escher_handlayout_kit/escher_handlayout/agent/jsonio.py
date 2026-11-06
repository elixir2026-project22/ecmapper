"""Extract exactly one JSON object from a model reply - strictly, without executing anything."""
import re
from ..declspec import SpecError, parse_json_strict

MAX_REPLY_CHARS = 400_000

def _first_object(text):
    start = text.find('{')
    while start != -1:
        depth, in_str, esc = 0, False, False
        for i in range(start, len(text)):
            c = text[i]
            if in_str:
                if esc: esc = False
                elif c == '\\': esc = True
                elif c == '"': in_str = False
            else:
                if c == '"': in_str = True
                elif c == '{': depth += 1
                elif c == '}':
                    depth -= 1
                    if depth == 0: return text[start:i + 1]
        start = text.find('{', start + 1)
    return None

def extract_json_object(text):
    """Return the parsed first top-level JSON object found in `text` (code fences and surrounding prose tolerated)."""
    if not isinstance(text, str) or not text.strip(): raise SpecError('empty reply')
    if len(text) > MAX_REPLY_CHARS: raise SpecError('reply too long')
    cands = [m.group(1) for m in re.finditer(r'```(?:json|JSON)?\s*(.*?)```', text, re.S)] + [text]
    for c in cands:
        obj = _first_object(c)
        if obj is not None:
            val = parse_json_strict(obj)
            if not isinstance(val, dict): raise SpecError('the reply must be a JSON object')
            return val
    raise SpecError('no JSON object found in the reply')
