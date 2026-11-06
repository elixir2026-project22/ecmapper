#!/usr/bin/env python3
"""Model-agnostic layout agent: a language model proposes layout spec fragments (JSON) in a chat; trusted code validates and applies them.

Examples
  # no API at all: copy/paste each prompt into any chat model (ChatGPT, Claude, Gemini, a local model ...)
  python agent.py --input map.json --out run1 --provider manual

  # OpenAI-compatible endpoint (OpenAI, OpenRouter, Ollama, vLLM, LM Studio ...). The key is read from an environment variable.
  export OPENAI_API_KEY=...
  python agent.py --input map.json --out run1 --provider openai --model gpt-4o
  python agent.py --input map.json --out run1 --provider openai --base-url http://localhost:11434/v1 --model llama3.1 --api-key-env ""

  # Anthropic Messages API
  export ANTHROPIC_API_KEY=...
  python agent.py --input map.json --out run1 --provider anthropic --model <model-name>

Use --resume to continue an interrupted run (state is saved in <out>/state.json).
The model never executes code: its replies are parsed as JSON only (see docs/AGENT.md for the security model).
"""
import argparse, os, sys, json
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from escher_handlayout.agent.session import Session, AgentFailure
from escher_handlayout.agent import llm as backends
from escher_handlayout._io import read_json

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('--input', required=True, help='draft Escher map (.json)')
ap.add_argument('--out', required=True, help='output directory (created)')
ap.add_argument('--provider', choices=['manual', 'openai', 'anthropic', 'scripted'], default='manual')
ap.add_argument('--model', help='model name (openai/anthropic providers)')
ap.add_argument('--base-url', help='endpoint base URL (openai-compatible: e.g. https://api.openai.com/v1, http://localhost:11434/v1)')
ap.add_argument('--api-key-env', help='NAME of the environment variable holding the API key (default OPENAI_API_KEY / ANTHROPIC_API_KEY; "" = no key)')
ap.add_argument('--temperature', type=float, default=0.2); ap.add_argument('--max-tokens', type=int, default=8000)
ap.add_argument('--json-mode', action='store_true', help="ask the endpoint for JSON-only output (response_format json_object), if it supports it")
ap.add_argument('--script', help='JSON file with a list of prepared replies (provider scripted; for tests)')
ap.add_argument('--curv', type=float, default=0.6); ap.add_argument('--curv-candidates', default='0.3,0.45,0.6')
ap.add_argument('--max-refine', type=int, default=4); ap.add_argument('--max-retries', type=int, default=3)
ap.add_argument('--max-turns', type=int, default=80); ap.add_argument('--chunk', type=int, default=40, help='max reactions per layout turn')
ap.add_argument('--resume', action='store_true')
ap.add_argument('--seed-spec', help='JSON spec (fragment) to start from, e.g. your own zones/plan/merges - or a complete spec; the model only fills in what is missing')
a = ap.parse_args()

if a.provider == 'manual': llm = backends.Manual(os.path.join(a.out, 'manual'))
elif a.provider == 'scripted': llm = backends.Scripted(read_json(a.script))
else:
    env = a.api_key_env if a.api_key_env is not None else ('OPENAI_API_KEY' if a.provider == 'openai' else 'ANTHROPIC_API_KEY')
    key = os.environ.get(env) if env else None
    local = (a.base_url or '').startswith(('http://localhost', 'http://127.0.0.1'))
    if env and not key and not local: sys.exit(f'environment variable {env} is not set (the key is never taken from the command line)')
    if not a.model: sys.exit('--model is required for this provider')
    if a.provider == 'openai':
        llm = backends.OpenAICompat(a.base_url or 'https://api.openai.com/v1', a.model, key, a.temperature, a.max_tokens, a.json_mode)
    else:
        llm = backends.Anthropic(a.model, key, a.base_url or 'https://api.anthropic.com', a.temperature, a.max_tokens)
cands = tuple(float(x) for x in a.curv_candidates.split(',') if x.strip())
sess = Session(a.input, a.out, llm, curv=a.curv, curv_candidates=cands, max_refine=a.max_refine, max_retries=a.max_retries,
               max_turns=a.max_turns, chunk=a.chunk, resume=a.resume,
               seed_spec=read_json(a.seed_spec) if a.seed_spec else None)
try:
    s = sess.run()
except (AgentFailure, backends.LLMError) as e:
    print('STOPPED:', e); sys.exit(3)
m = s['metrics']
print(f"\nDONE. map: {s['outputs']['map']}\n  png: {s['outputs']['png']}\n  turns: {s['turns']}  curvature: {s['curvature']}  targets met: {s['targets_met']}")
print(f"  crossings={m['line_crossings']} labels_on_lines={m['labels_on_lines']} labels_on_nodes={m['labels_on_nodes']} overlaps={m['label_label_overlaps']}  validation errors: {len(s['validation_errors'])}")
