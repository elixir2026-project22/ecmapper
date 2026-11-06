"""Chat session that drives a language model through plan -> layout -> completion -> refinement -> finalisation.

The model only returns JSON spec fragments. Every fragment is parsed strictly, schema-checked, merged into the spec and
checked against the map before it is accepted; errors go back to the model as the next chat message. All map data shown to the
model is sanitised. Nothing from a model reply is executed, imported or used as a path.
"""
import os, json, hashlib, re, copy, shutil, time
from ..engine import Layout
from ..declspec import SpecError, compile_spec, merge_fragment, validate_schema
from ..inspect_map import analyse
from .. import pipeline
from .jsonio import extract_json_object
from .llm import LLMError
from .._io import read_text, read_json, write_json, write_text

PROMPT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'prompts')


class AgentFailure(Exception):
    pass


def read_prompt(name):
    return read_text(os.path.join(PROMPT_DIR, name))

def fill(t, **kw):
    for k, v in kw.items(): t = t.replace('{{' + k + '}}', str(v))
    return t

def safe(s, n=40):
    """Sanitise text taken from the map before it is shown to a model (no control chars, no markup, short)."""
    s = re.sub(r'[^A-Za-z0-9 _.,:;()+\-/\']', '?', str(s)).strip()
    return s[:n]

def system_prompt():
    return fill(read_prompt('system.md'), FORMAT_REFERENCE=read_prompt('format_reference.md'), RECIPES=read_prompt('recipes.md'))

def score(rep):
    m = rep['metrics']
    return (1000 * len(rep['validation_errors']) + m['line_crossings'] + 3 * m['labels_on_lines'] + 6 * m['labels_on_nodes']
            + 4 * m['label_label_overlaps'] + 10 * m['primary_nodes_closer_than_70'])

def targets_met(rep):
    m = rep['metrics']; n = rep['n_reactions']
    return (not rep['validation_errors'] and m['labels_on_nodes'] == 0 and m['label_label_overlaps'] <= 2
            and m['primary_nodes_closer_than_70'] == 0 and m['labels_on_lines'] <= max(6, 0.07 * n) and m['line_crossings'] <= 0.35 * n)

def zone_name(title):
    s = re.sub(r"[^A-Za-z0-9 .,()/+:'_\-]", ' ', title).strip()
    s = re.sub(r'\s+', ' ', s)[:80] or 'Module'
    return s[0].upper() + s[1:] if s[0].isalpha() else 'M ' + s


class Session:
    def __init__(self, input_path, out_dir, llm, curv=0.6, curv_candidates=(0.3, 0.45, 0.6), max_refine=4, max_retries=3,
                 max_turns=80, chunk=40, resume=False, log=print, seed_spec=None):
        self.input = os.path.abspath(input_path); self.out = os.path.abspath(out_dir); os.makedirs(self.out, exist_ok=True)
        self.llm = llm; self.curv = curv; self.cands = tuple(curv_candidates); self.max_refine = max_refine
        self.max_retries = max_retries; self.max_turns = max_turns; self.chunk = chunk; self.log = log
        self.system = system_prompt()
        self.name = os.path.splitext(os.path.basename(self.input))[0] + '_handlayout'
        h = hashlib.sha256(read_text(self.input, 'rb')).hexdigest()
        self.state_path = os.path.join(self.out, 'state.json')
        if resume and os.path.exists(self.state_path):
            self.state = read_json(self.state_path)
            if self.state.get('input_sha') != h: raise AgentFailure('state.json belongs to a different input map')
        else:
            self.state = dict(input_sha=h, step='plan', spec={'spec_version': 1}, tasks=[], task_i=0, refine_i=0, stalls=0,
                              best=None, turns=0, scores=[], notes=[])
        os.makedirs(os.path.join(self.out, 'transcript'), exist_ok=True)
        if seed_spec and not (resume and os.path.exists(self.state_path)): self._seed(seed_spec)

    # ------------------------------------------------------------------ plumbing
    def save(self):
        tmp = self.state_path + '.tmp'; write_json(tmp, self.state); os.replace(tmp, self.state_path)

    @property
    def spec(self): return self.state['spec']

    def fresh(self): return Layout(self.input)

    def _turn(self, kind, messages):
        if self.state['turns'] >= self.max_turns: raise AgentFailure(f'turn budget of {self.max_turns} exhausted')
        self.state['turns'] += 1; n = self.state['turns']
        self.log(f'[{kind}] turn {n} ...')
        reply = self.llm.chat(messages, turn=n)
        write_json(os.path.join(self.out, 'transcript', f'{n:03d}_{kind}.json'), dict(kind=kind, turn=n, messages=messages[-3:] if n > 1 else messages, reply=reply))
        return reply

    def ask(self, kind, user_text, validator):
        """One chat exchange with automatic retries: returns (fragment, validator_result)."""
        messages = [{'role': 'system', 'content': self.system}, {'role': 'user', 'content': user_text}]
        last = ''
        for attempt in range(self.max_retries + 1):
            text = self._turn(kind if attempt == 0 else kind + '_retry', messages)
            try:
                frag = extract_json_object(text)
                errs = validate_schema(frag)
                if errs: raise SpecError('schema problems:\n  ' + '\n  '.join(errs[:20]))
                return frag, validator(frag)
            except SpecError as e:
                last = str(e); self.log(f'    rejected: {last.splitlines()[0][:150]}')
                messages += [{'role': 'assistant', 'content': text[:6000]}, {'role': 'user', 'content': fill(read_prompt('fix.md'), ERRORS=last[:3000])}]
        raise AgentFailure(f'the model did not produce an acceptable {kind} reply after {self.max_retries + 1} attempts. Last problem: {last[:500]}')


    def check_plan(self, new, allr):
        """Problems of the plan section (list of strings): complete, consistent partition of the reaction ids."""
        errs, seen = [], {}
        names = {z['name'] for z in new.get('zones', [])}
        for z, e in new.get('plan', {}).items():
            if z not in names: errs.append(f"plan zone '{z}' is not in zones")
            if 'origin' not in e or 'size' not in e or 'reactions' not in e: errs.append(f"plan['{z}'] needs origin, size and reactions")
            for r in e.get('reactions', []):
                if r not in allr: errs.append(f"plan['{z}']: reaction id {r} does not exist")
                elif r in seen: errs.append(f"reaction {r} is listed in zones '{seen[r]}' and '{z}'")
                seen[r] = z
        miss = sorted(allr - set(seen), key=int)
        if miss: errs.append(f'{len(miss)} reactions are in no zone: ' + ', '.join(miss[:40]))
        for z in names:
            if z not in new.get('plan', {}): errs.append(f"zone '{z}' has no plan entry")
        return errs

    def make_tasks(self, new):
        tasks = []
        for z, e in new['plan'].items():
            rids = e['reactions']
            for i in range(0, len(rids), self.chunk):
                tasks.append(dict(zone=z, block=z if i == 0 else f'{z} ({i // self.chunk + 1})', reactions=rids[i:i + self.chunk],
                                  first=(i == 0), part=i // self.chunk + 1, parts=(len(rids) + self.chunk - 1) // self.chunk))
        return tasks

    def _seed(self, seed):
        errs = validate_schema(seed)
        if errs: raise AgentFailure('seed spec violates the schema: ' + '; '.join(errs[:10]))
        new = merge_fragment({'spec_version': 1}, seed)
        L = self.fresh(); P, W = compile_spec(new, L, partial=True)
        if P: raise AgentFailure('seed spec problems: ' + '; '.join(P[:10]))
        self.state['spec'] = new
        if new.get('plan'):
            errs = self.check_plan(new, set(L.R))
            if errs: raise AgentFailure('seed plan problems: ' + '; '.join(errs[:10]))
            self.state['tasks'] = self.make_tasks(new); self.state['step'] = 'layout'; self.state['task_i'] = 0
        self.save()

    # ------------------------------------------------------------------ data shown to the model
    def view(self, spec=None):
        """(Layout after a partial compile of the spec, alias->id map, id->alias map)."""
        L = self.fresh(); compile_spec(spec or self.spec, L, partial=True)
        return L, L.spec_alias, {v: k for k, v in L.spec_alias.items()}

    def rx_row(self, L, rid, rev_alias):
        r = L.R[rid]; I = L.info[rid]
        def f(n):
            nm = safe(L.N[n]['bigg_id'], 24)
            return f"{nm} [{rev_alias.get(n, n)}*]" if L.N[n]['node_is_primary'] else nm
        return (f"{rid} | {safe(r['bigg_id'], 24)} ({safe(r.get('name', ''), 36)}) | {'<->' if r['reversibility'] else '->'} | "
                + ' + '.join(f(n) for _, n in I['subs']) + ' => ' + ' + '.join(f(n) for _, n in I['prods']))

    def propose(self, analysis):
        comps = analysis['components']; zones = {}; order = []
        for r in analysis['reactions']:
            c = comps[r['component']]; title = c['suggested_subsystem'] or f"Module {c['index'] + 1}"
            if title not in zones: zones[title] = []; order.append(title)
            zones[title].append(r['id'])
        out = dict(zones=[dict(name=zone_name(t), title=zone_name(t), label=t if any(x['text'] == t for x in analysis['title_labels']) else None) for t in order],
                   reactions_by_zone={zone_name(t): zones[t] for t in order})
        return out

    # ------------------------------------------------------------------ steps
    def step_plan(self):
        analysis = analyse(self.input)
        write_json(os.path.join(self.out, 'analysis.json'), analysis)
        if analysis['unsupported_reactions']:
            raise AgentFailure('unsupported reaction topology (no multimarkers): ' + ', '.join(r['bigg_id'] for r in analysis['unsupported_reactions'][:20]))
        prop = self.propose(analysis); self.proposal = prop
        overview = [f"{analysis['n_reactions']} reactions, {analysis['n_nodes']} nodes, {len(analysis['components'])} connected components."]
        for z in prop['zones']:
            rids = prop['reactions_by_zone'][z['name']]
            names = [safe(r['bigg_id'], 18) for r in analysis['reactions'] if r['id'] in set(rids)][:40]
            overview.append(f"- {z['name']}: {len(rids)} reactions ({', '.join(names)}{' ...' if len(rids) > 40 else ''})")
        rows = ['| metabolite | node ids | reactions (first) | hub |', '|---|---|---|---|']
        for c in analysis['merge_candidates'][:60]:
            rows.append(f"| {safe(c['bigg_id'], 20)} | {', '.join(c['nodes'])} | {safe(', '.join(c['reactions'][:6]), 60)} | {'hub' if c['in_default_hub_list'] else ''} |")
        user = fill(read_prompt('plan.md'), OVERVIEW='\n'.join(overview), PROPOSAL=json.dumps(prop, indent=1)[:12000], MERGE_TABLE='\n'.join(rows))
        allr = {r['id'] for r in analysis['reactions']}
        def validate(frag):
            new = merge_fragment(self.spec, frag)
            if not new.get('zones') or not new.get('plan'): raise SpecError('the plan must contain "zones" and "plan"')
            P, W = compile_spec(new, self.fresh(), partial=True)
            if P: raise SpecError('\n'.join(P[:25]))
            errs = self.check_plan(new, allr)
            if errs: raise SpecError('\n'.join(errs[:30]))
            return new
        frag, new = self.ask('plan', user, validate)
        st = self.state; st['spec'] = new
        st['tasks'] = self.make_tasks(new)
        if frag.get('notes'): st['notes'].append(frag['notes'])
        st['step'] = 'layout'; st['task_i'] = 0; self.save()

    def step_layout(self):
        st = self.state
        while st['task_i'] < len(st['tasks']):
            t = st['tasks'][st['task_i']]
            L, alias, rev = self.view()
            if all(r in L.spec_done for r in t['reactions']): st['task_i'] += 1; continue      # already placed (seed spec)
            plan = self.spec['plan'][t['zone']]
            used = []
            for rid in t['reactions']:
                for _, n in L.info[rid]['subs'] + L.info[rid]['prods']:
                    if L.N[n]['node_is_primary'] and n not in used: used.append(n)
            to_place = [n for n in used if n not in L.pos]; placed = [n for n in used if n in L.pos]
            fmtn = lambda n: f"[{rev.get(n, n)}] {safe(L.N[n]['bigg_id'], 24)}"
            tp = '\n'.join(fmtn(n) for n in to_place) or '(none - all nodes of these reactions are already placed)'
            pl = '\n'.join(f"{fmtn(n)} at ({L.pos[n][0]:.0f}, {L.pos[n][1]:.0f})" for n in placed) or '(none)'
            zb = {}
            for nid, bn in L.spec_placed.items():
                z = next((b.get('zone') for b in self.spec.get('layouts', []) if b['block'] == bn), None)
                zb.setdefault(z, []).append(L.pos[nid])
            occ = [f"zone {z}: x {min(p[0] for p in v):.0f}..{max(p[0] for p in v):.0f}, y {min(p[1] for p in v):.0f}..{max(p[1] for p in v):.0f}" for z, v in zb.items()]
            for z, e in self.spec['plan'].items():
                if z not in zb and z != t['zone']: occ.append(f"zone {z} (not yet laid out): planned x {e['origin'][0]}..{e['origin'][0] + e['size'][0]}, y {e['origin'][1]}..{e['origin'][1] + e['size'][1]}")
            first_chunk = st['tasks'][[i for i, x in enumerate(st['tasks']) if x['zone'] == t['zone']][0]]
            user = fill(read_prompt('layout.md'), ZONE=t['zone'], BLOCK=t['block'], ORIGIN=json.dumps(plan['origin']), SIZE=json.dumps(plan['size']),
                        CHUNK_INFO=f"part {t['part']} of {t['parts']}, {len(t['reactions'])} reactions", OCCUPIED='; '.join(occ) or 'nothing yet',
                        TITLE_HINT=(f', "titles": {{"{t["zone"]}": [x, y]}}' if t['first'] else ''),
                        TO_PLACE=tp, PLACED=pl, REACTIONS='\n'.join(self.rx_row(L, r, rev) for r in t['reactions']))
            want = set(t['reactions'])
            def validate(frag, t=t, want=want):
                new = merge_fragment(self.spec, frag)
                L2 = self.fresh(); P, W = compile_spec(new, L2, partial=True)
                if P: raise SpecError('\n'.join(P[:25]))
                miss = [r for r in want if r not in L2.spec_done]
                if miss:
                    why = []
                    for r in miss[:20]:
                        lack = [L2.N[n]['bigg_id'] + '#' + n for _, n in L2.info[r]['subs'] + L2.info[r]['prods'] if L2.N[n]['node_is_primary'] and n not in L2.pos]
                        why.append(f"{r} ({L2.R[r]['bigg_id']})" + (f" lacks coordinates for {', '.join(lack)}" if lack else ' is not described or has unresolved references'))
                    raise SpecError(f'{len(miss)} of the requested reactions are not placed: ' + '; '.join(why) + '. Send the complete fragment again.')
                return new
            frag, new = self.ask('layout', user, validate)
            st['spec'] = new
            if frag.get('notes'): st['notes'].append(frag['notes'])
            st['task_i'] += 1; self.save()
        st['step'] = 'complete'; self.save()

    def step_complete(self):
        """Every reaction must be placed before a full build is possible; ask for what is missing."""
        st = self.state
        def validate(frag):
            new = merge_fragment(self.spec, frag)
            compile_spec(new, self.fresh())           # full check (raises SpecError)
            return new
        try:
            compile_spec(self.spec, self.fresh()); st['step'] = 'refine'; self.save(); return
        except SpecError as e:
            problem = str(e)
        L, alias, rev = self.view()
        miss = [r for r in L.R if r not in L.spec_done]
        user = ("## TASK: COMPLETE the layout\n\nThe full check found these problems:\n\n" + problem[:3500] +
                "\n\n### DATA: the reactions that are still unplaced\n" + '\n'.join(self.rx_row(L, r, rev) for r in miss[:80]) +
                "\n\nReply with ONE JSON fragment that fixes them (blocks/nodes/reactions to add or correct; reuse the existing block names and origins).")
        frag, new = self.ask('complete', user, validate)
        st['spec'] = new; st['step'] = 'refine'; self.save()

    def build(self, spec, curv, tag, png=False, name=None):
        d = os.path.join(self.out, 'iterations', tag)
        return pipeline.build(self.input, spec, d, curv=curv, do_pack=True, png=png, name=name or 'iter')

    def offenders(self, rep):
        m = rep['metrics']; names = []
        for (a, b), n in m['worst_crossing_pairs']: names += [a, b] * n
        for s in m['labels_on_lines_list'] + m['labels_on_nodes_list']:
            if s.startswith('reaction '): names.append(s[9:])
        for a, b in m['label_label_overlaps_list']:
            for s in (a, b):
                if s.startswith('reaction '): names.append(s[9:])
        counts = {}
        for n in names: counts[n] = counts.get(n, 0) + 1
        return [n for n, _ in sorted(counts.items(), key=lambda kv: -kv[1])][:12]

    def entries(self, names, spec):
        L, _, _ = self.view(spec); out = []
        bigg = {}
        for rid, r in L.R.items(): bigg.setdefault(r['bigg_id'], []).append(rid)
        want = {rid for n in names for rid in bigg.get(n, [])}
        for b in spec.get('layouts', []):
            for r in b.get('reactions', []):
                if r['id'] in want:
                    prim = {n: [round(L.pos[n][0]), round(L.pos[n][1])] for _, n in L.info[r['id']]['subs'] + L.info[r['id']]['prods'] if L.N[n]['node_is_primary'] and n in L.pos}
                    out.append(dict(block=b['block'], origin=b.get('origin', [0, 0]), entry=r, bigg_id=L.R[r['id']]['bigg_id'], primary_nodes_global=prim))
        return out

    def step_refine(self):
        st = self.state
        if st['best'] is None:
            rep = self.build(self.spec, self.curv, 'it_00'); st['best'] = dict(score=score(rep), spec=copy.deepcopy(self.spec), rep=rep)
            st['scores'].append(('initial', st['best']['score'])); self.log(f"[build] initial score {st['best']['score']}  crossings={rep['metrics']['line_crossings']}"); self.save()
        while st['refine_i'] < self.max_refine and st['stalls'] < 2:
            rep = st['best']['rep']
            if targets_met(rep): self.log('[refine] targets met'); break
            names = self.offenders(rep)
            if not names: break
            body = pipeline.report_md(rep)
            user = fill(read_prompt('refine.md'), ROUND=st['refine_i'] + 1, MAX_ROUND=self.max_refine, REPORT=body[:3500],
                        ENTRIES=json.dumps(self.entries(names, st['best']['spec']), indent=1)[:9000])
            def validate(frag):
                new = merge_fragment(st['best']['spec'], frag)
                compile_spec(new, self.fresh())
                return new
            frag, new = self.ask('refine', user, validate)
            st['refine_i'] += 1
            rep2 = self.build(new, self.curv, f"it_{st['refine_i']:02d}"); sc = score(rep2)
            st['scores'].append((f"refine {st['refine_i']}", sc)); self.log(f"[refine] round {st['refine_i']}: score {sc} (best {st['best']['score']})")
            if sc < st['best']['score']:
                st['best'] = dict(score=sc, spec=copy.deepcopy(new), rep=rep2); st['stalls'] = 0
            else:
                st['stalls'] += 1
            self.save()
        st['spec'] = copy.deepcopy(st['best']['spec']); st['step'] = 'finalize'; self.save()

    def step_finalize(self):
        st = self.state; best = None
        for c in sorted(set(self.cands) | {self.curv}):
            rep = self.build(st['spec'], c, f'curv_{c:.2f}'); sc = score(rep)
            st['scores'].append((f'curvature {c:.2f}', sc)); self.log(f'[finalize] curvature {c}: score {sc}')
            key = (sc, abs(c - self.curv))
            if best is None or key < best[0]: best = (key, c)
        c = best[1]
        rep = pipeline.build(self.input, st['spec'], self.out, curv=c, do_pack=True, png=True, name=self.name)
        write_json(os.path.join(self.out, self.name + '.layout.json'), st['spec'])
        summary = dict(model=self.llm.name, turns=st['turns'], curvature=c, scores=st['scores'], notes=st['notes'],
                       validation_errors=rep['validation_errors'], metrics={k: v for k, v in rep['metrics'].items() if not k.endswith('list')},
                       area_before=rep['area_before_packing_Munits2'], area_after=rep['area_after_packing_Munits2'], canvas=rep['canvas'],
                       targets_met=targets_met(rep), outputs=dict(map=rep['output_json'], png=rep.get('output_png'), spec=os.path.join(self.out, self.name + '.layout.json')))
        write_json(os.path.join(self.out, 'run_summary.json'), summary)
        st['step'] = 'done'; st['summary'] = summary; self.save()

    # ------------------------------------------------------------------ driver
    def run(self):
        try:
            while self.state['step'] != 'done': getattr(self, 'step_' + self.state['step'])()
        except (AgentFailure, LLMError) as e:
            self.save()
            write_text(os.path.join(self.out, 'FAILED.md'), f"# Run stopped\n\n{e}\n\nThe session state was saved; rerun with --resume to continue.\n")
            raise
        return self.state['summary']
