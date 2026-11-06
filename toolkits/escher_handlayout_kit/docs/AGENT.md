# The chat-based layout agent (`agent.py`)

A small program that lets **any language model** lay out an Escher draft map in a chat - while the model is **never able to run code**.

```
  you / config ──► agent.py (trusted code) ◄──────── text only ────────► language model (untrusted)
                     │  1. reads the map, builds prompts (data is sanitised)
                     │  2. sends one chat turn, gets back ONE JSON object
                     │  3. parses it strictly, checks it against a schema AND against the map
                     │  4. merges it into the layout spec (a JSON file) - or sends the error back to the model
                     └─ 5. trusted code (escher_handlayout) builds geometry, curves, labels, packing, QA metrics
```

The model's only job is the *creative* part: where each subsystem goes, where the nodes are, which side the cofactors hang on.
Everything mechanical (curves, label placement, packing, validation) is deterministic code that you control.

## 1. Three ways to run it

### A. Copy/paste mode - no API, no keys, any chat window (recommended to start)

```bash
pip install -r requirements.txt
python agent.py --input my_draft_map.json --out run1 --provider manual
```
What happens: the program writes `run1/manual/turn_001_prompt.md`. Open it, copy **everything** into a chat with any model
(ChatGPT, Claude, Gemini, a local model ...), copy the model's answer, paste it into the terminal and finish with a line
containing only `END` (or save it as `run1/manual/turn_001_reply.txt` and press Enter). The program checks the answer,
sometimes sends you a new prompt with an error message to paste back (that is normal), and continues. A 90-reaction map takes
about 10-20 turns. Stop any time; continue later with `--resume`.

### B. Automatic mode - OpenAI-compatible API (OpenAI, OpenRouter, Ollama, vLLM, LM Studio, llama.cpp ...)

```bash
export OPENAI_API_KEY=...                       # the key is read from an environment variable, never from the command line
python agent.py --input map.json --out run1 --provider openai --model <model-name>
# local model (no key):
python agent.py --input map.json --out run1 --provider openai --base-url http://localhost:11434/v1 --model <name> --api-key-env ""
```
### C. Automatic mode - Anthropic API
```bash
export ANTHROPIC_API_KEY=...
python agent.py --input map.json --out run1 --provider anthropic --model <model-name>
```
Useful options: `--json-mode` (ask the endpoint for JSON-only output), `--temperature 0.2`, `--max-retries 3`, `--max-refine 4`,
`--max-turns 80`, `--chunk 40` (reactions per layout turn; use 20 for weaker models), `--curv 0.6`, `--resume`.

## 2. What the conversation looks like

| turn | harness sends | model returns (one JSON object) |
|---|---|---|
| 1 plan | overview of the map, proposed subsystems, hub-merge candidates | `zones`, `plan` (region + reaction list per zone), `merges`, `splits` |
| 2..n layout | one zone (or a chunk of 40 reactions): nodes to place, already placed hubs, free space, reaction table | one `layouts` block (+ the zone `titles`) |
| complete | only if something is still unplaced: the list of problems | the missing entries |
| refine (<=4) | quality report with named offenders (crossings, labels on lines ...) and their current entries | a *small correction*, or `{"notes": "no change"}` |
| finalize | - (no model) | the program tries curvature 0.3 / 0.45 / 0.6, keeps the best, writes the files |

Every turn is a fresh, short chat (system prompt about 2,100 tokens + 1,000-1,400 tokens of data). Measured on the bundled
91-reaction example: 11 turns, about 38k input tokens and 3k output tokens without refinement. Retries after a rejected reply
append the error message to that turn's chat.

Outputs in `--out`: `<map>_handlayout.json` (the Escher map), `..._overview.png`, `qa_report.md/json`, `<map>_handlayout.layout.json`
(the final spec - reusable with `python run.py build ...`), `run_summary.json`, `transcript/` (every prompt and reply), `state.json`.

## 3. Security model

| risk | what prevents it |
|---|---|
| **The model executes code** | It cannot: it has no tools, and its reply is only ever parsed by a strict JSON parser. There is no `eval`, `exec`, import or file access based on model text; a test scans the agent code for such calls. Python specs (`.py`) are disabled by default (`--unsafe-python-spec` exists only for specs *you* wrote). |
| **Malicious/odd JSON** | Strict schema (unknown keys rejected, enums, numeric ranges +-20000, id/name patterns without control characters or angle brackets), duplicate keys and NaN/Infinity rejected, size caps (reply 400k chars, spec 1.5 MB, counts per section). |
| **Nonsense or hostile content in the map** (reaction names ...) | Only ids and sanitised short names are shown to the model, inside clearly marked DATA sections; the system prompt tells the model to treat them as data. Even a model that follows an injected instruction can only move nodes or merge nodes of the *same* metabolite. |
| **Biochemistry altered** | Not expressible in the format; the final validator also compares reactions, stoichiometry, reversibility and genes of output and input. |
| **Resource exhaustion** | Turn budget (`--max-turns`), retry limit, bounded geometry (coordinates, counts), packing fallback for oversized pieces. |
| **Secrets** | API key only from an environment variable (never a CLI argument, never written to disk); transcripts contain prompts and replies only. |
| **Paths** | The model never supplies a file name or path. |

Residual risks / things you should know:
- **Privacy:** the map content (ids, names, structure) is sent to the model provider in the prompts. For unpublished models use a local model or manual mode with a provider you trust.
- A model can still produce an *ugly* layout. The numeric QA metrics catch crossings, labels on lines etc., not aesthetics.
- For the strictest setup run `agent.py` in a container with network access only to your model endpoint (the program itself needs no other network).

## 4. Choosing and tuning models

- Any model that can follow a JSON format works in principle; **spatial reasoning is the hard part**. Stronger models place nodes
  more sensibly. Weak models: `--chunk 20`, `--max-retries 5`, `--json-mode` where supported, lower temperature.
- The prompt files in `escher_handlayout/agent/prompts/` are plain Markdown - edit them to tune wording for a particular model.
- `plan` decides the global arrangement (the model sees only reaction ids/names, not a picture). If you want a specific
  arrangement, write the `plan` yourself: save a spec with `zones` + `plan` + `merges`, and start the session with it
  (copy it into `<out>/state.json` as `"spec"` with `"step": "layout"` - see `Session` in `agent/session.py`) - or simply give the
  instruction in manual mode by adding a sentence to the first prompt.

## 5. Troubleshooting

| symptom | meaning / remedy |
|---|---|
| "rejected: schema problems ..." | the model used a wrong key/value; the message is sent back automatically; if it repeats, lower temperature or use a stronger model |
| "node X was merged into Y" | after merges the old id no longer exists; the model must use the `keep` id |
| "N of the requested reactions are not placed" | the model forgot reactions or nodes; it is asked again |
| "the model did not produce an acceptable reply" | `FAILED.md` and `state.json` are written; fix the problem (e.g. other model) and `--resume` |
| QA metrics stay high after refinement | normal for dense maps; the best build is kept; inspect `qa_report.md` and edit the final `.layout.json` by hand |

## 6. Status and limits (honest)

- Tested offline: a scripted model replays the example layout turn by turn and reproduces the expected output exactly; invalid
  replies, retries, budget/resume, copy/paste mode, malformed JSON, schema and semantic rejections, random-spec fuzzing and the
  HTTP clients (against a local fake server) are covered by `python tests/run_tests.py`.
- **Not yet tested with real language models.** Prompt wording will probably need tuning for specific models; the transcripts make
  that easy to see.
- Works for maps with the Escher *draft-map topology* only (see README).
- The OpenAI/Anthropic clients use only the standard library and cover the plain chat-completion calls (no streaming, no tools).
