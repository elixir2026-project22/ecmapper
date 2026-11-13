"""Score and targets of a build, the history of the builds of a job, and what to do next.

The build prints all of this, so that whoever iterates on a spec (a person or a model) neither has to compute the
score nor to decide alone when to stop (PROMPT.md step E).
"""
import hashlib, json, os
from ._io import read_json, write_json, write_text

#: weight of each term of the score, in the order they are printed
WEIGHTS = (('validation errors', 1000), ('crossings', 1), ('labels on lines', 3), ('labels on nodes', 6), ('label overlaps', 4),
           ('close nodes', 10), ('unconnected subsystems', 40), ('detached pieces', 20), ('placement penalty', 1))
#: what to change in the spec for a missed target (docs/GUIDE.md section 8)
REMEDY = {
    'validation errors': 'read the ERROR lines of qa_report.md: usually a missing placement or a wrong merge',
    'crossings': 'for the worst pairs: move one node, flip a `side`, change a lane `shift`, or split a node instead of drawing a long line; '
                 'if many pairs involve two subsystems, move or mirror one block (`origin`)',
    'labels on lines': 'widen `dist` / `step` / the lane spacing of the named reactions',
    'labels on nodes': 'move the node or give its neighbours more room',
    'label overlaps': 'spread the two reactions apart (150 units or more) or change a cofactor `side`',
    'close nodes': 'two primary nodes are nearly on top of each other: correct their coordinates (merge them if they are one metabolite)',
    'unconnected subsystems': 'connect ONE subsystem per build, the first one named: move its block (`origin`) next to the subsystem it shares the '
                              'metabolite with, on the side analysis/arrangement.md says, and patch the named merge (the list of dups is complete: '
                              'send it as it is). Then build. Do not tidy the inside of a block before all subsystems are connected',
    'detached pieces': 'merge the shared metabolite into the main map (`merges`), or keep the copy on purpose and say why in DECISIONS.md',
    'placement penalty': 'move the block (`origin`) where the direction says; a medium direction may stay violated if following it costs crossings',
}


def terms(rep):
    """The numbers the score is made of: [(name, value, weight)]."""
    m = rep['metrics']
    values = {'validation errors': len(rep['validation_errors']), 'crossings': m['line_crossings'], 'labels on lines': m['labels_on_lines'],
              'labels on nodes': m['labels_on_nodes'], 'label overlaps': m['label_label_overlaps'],
              'close nodes': m['primary_nodes_closer_than_70'], 'unconnected subsystems': len(rep.get('unconnected') or []),
              'detached pieces': len(rep.get('detached') or []),
              'placement penalty': (rep.get('placement') or {}).get('penalty', 0)}
    return [(name, values[name], w) for name, w in WEIGHTS]


def score(rep):
    return sum(v * w for _, v, w in terms(rep))


def score_line(rep):
    parts = [f'{v} {name}' if w == 1 else f'{w}*{v} {name}' for name, v, w in terms(rep) if v]
    return f"score {score(rep)}" + (' = ' + ' + '.join(parts) if parts else '') + '  (lower is better)'


def targets(rep):
    """The targets of PROMPT.md: [dict(name, value, limit, ok)]. 'detached pieces' has ok None: a piece may be a deliberate copy."""
    m = rep['metrics']; n = rep['n_reactions']; pl = rep.get('placement') or {}
    strong = sum(1 for s in pl.get('violated_strengths', []) if s in ('hard', 'strong'))
    rows = [('validation errors', len(rep['validation_errors']), 0), ('crossings', m['line_crossings'], int(0.1 * n)),
            ('labels on lines', m['labels_on_lines'], int(max(6, 0.07 * n))), ('labels on nodes', m['labels_on_nodes'], 0),
            ('label overlaps', m['label_label_overlaps'], 2), ('close nodes', m['primary_nodes_closer_than_70'], 0),
            ('unconnected subsystems', len(rep.get('unconnected') or []), 0), ('violated hard/strong directions', strong, 0)]
    out = [dict(name=k, value=v, limit=lim, ok=v <= lim) for k, v, lim in rows]
    out.append(dict(name='detached pieces', value=len(rep.get('detached') or []), limit=0, ok=None if rep.get('detached') else True))
    return out


def targets_met(rep):
    return all(t['ok'] is not False for t in targets(rep))


def offenders(rep, name):
    """The named offenders of a target, as text ('' if the report names none)."""
    m = rep['metrics']
    if name == 'crossings': return '; '.join(f'{a} x {b} ({n})' for (a, b), n in m['worst_crossing_pairs'][:8])
    if name == 'labels on lines': return ', '.join(m['labels_on_lines_list'][:10])
    if name == 'labels on nodes': return ', '.join(m['labels_on_nodes_list'][:10])
    if name == 'label overlaps': return '; '.join(f'{a} / {b}' for a, b in m['label_label_overlaps_list'][:6])
    if name == 'unconnected subsystems':
        def one(u):
            s = u['shares'][0]; more = ', '.join(x['metabolite'] for x in u['shares'][1:4])
            return (f"{u['zone']} (shares {s['metabolite']} with {s['zone']}: " + '{"keep": "%s", "dups": [%s]}' % (s['keep'], ', '.join('"%s"' % d for d in s['dups']))
                    + (f'; also shares {more}' if more else '') + ')')
        return '; '.join(one(u) for u in rep['unconnected'][:4]) + (f" ... and {len(rep['unconnected']) - 4} more" if len(rep['unconnected']) > 4 else '')
    if name == 'detached pieces': return '; '.join(f"{', '.join(d['reactions'][:4])} (shares {', '.join(d['shares'][:4])})" for d in rep['detached'][:6])
    if name in ('violated hard/strong directions', 'placement penalty'): return '; '.join((rep.get('placement') or {}).get('violations', [])[:5])
    if name == 'validation errors': return '; '.join(rep['validation_errors'][:3])
    return ''


def target_lines(rep):
    out = []
    for t in targets(rep):
        mark = {True: 'ok  ', False: 'MISS', None: 'look'}[t['ok']]
        out.append(f"  {mark}  {t['name']}: {t['value']} (target: {'0' if t['limit'] == 0 else 'at most ' + str(t['limit'])})")
    return out


def first_fix(rep):
    """The missed target that costs most, its offenders and the remedy (None if nothing is missed)."""
    weight = dict(WEIGHTS); weight['violated hard/strong directions'] = 5
    missed = [t for t in targets(rep) if t['ok'] is not True]
    if not missed: return None
    t = max(missed, key=lambda t: (t['ok'] is False, t['name'] == 'unconnected subsystems', t['value'] * weight[t['name']]))     # connect first, tidy afterwards
    key = 'placement penalty' if t['name'].startswith('violated') else t['name']
    who = offenders(rep, t['name'])
    return f"{t['name']} ({t['value']})" + (f': {who}' if who else '') + f'. Remedy: {REMEDY[key]}.'


# ---------------------------------------------------------------------------------------------- history of a job
def spec_sha(spec_path):
    try:
        with open(spec_path, 'rb') as f: return hashlib.sha256(f.read().replace(b'\r\n', b'\n')).hexdigest()[:16]
    except OSError:
        return None


def record(history_path, rep, out_dir, spec_path):
    """Add this build to the history file (a rebuild into the same folder replaces its entry); returns the history."""
    hist = read_json(history_path) if os.path.exists(history_path) else []
    out = os.path.relpath(out_dir, os.path.dirname(os.path.abspath(history_path))).replace(os.sep, '/')
    hist = [h for h in hist if h['out'] != out]
    hist.append(dict(out=out, spec=spec_sha(spec_path), curvature=rep['curvature'], score=score(rep), targets_met=targets_met(rep),
                     terms={name: v for name, v, _ in terms(rep)}, map=os.path.basename(rep['output_json'])))
    write_json(history_path, hist)
    names = [name for name, _ in WEIGHTS]
    rows = ['| build | curvature | ' + ' | '.join(names) + ' | score | targets |', '|---|---|' + '---|' * (len(names) + 2)]
    rows += [f"| {h['out']} | {h['curvature']} | " + ' | '.join(str(h['terms'].get(n, 0)) for n in names) + f" | {h['score']} | {'met' if h['targets_met'] else 'missed'} |" for h in hist]
    write_text(os.path.splitext(history_path)[0] + '.md', '# Builds\n\n' + '\n'.join(rows) + '\n')
    return hist


def verdict(hist, rep, max_specs=16):
    """What to do after the newest build of `hist`: (state, text). States: finished, limit, stalled, go_on."""
    met = [h for h in hist if h['targets_met']]
    best = min(met or hist, key=lambda h: (h['score'], hist.index(h)))
    where =f"{best['out']}/{best['map']} (score {best['score']}, curvature {best['curvature']})"
    versions = []                                    # the first build of every version of the spec, in order
    for h in hist:
        if h['spec'] not in [v['spec'] for v in versions]: versions.append(h)
    if best['targets_met']:
        return 'finished', f'TARGETS MET by {where}. Finish: try the other curvatures once (--curv 0.3 and 0.45), then validate the best map and hand it in.'
    if len(versions) > max_specs:
        return 'limit', f'STOP: {len(versions) - 1} changed versions of the spec were built. Validate the best map and hand it in: {where}.'
    stalled = len(versions) >= 3 and all(v['score'] >= min(x['score'] for x in versions[:versions.index(v)]) for v in versions[-2:])
    fix = first_fix(rep)
    if stalled:
        return 'stalled', (f'NO IMPROVEMENT in the last two versions of the spec (best: {where}). Go back to the spec of the best build. If you have not '
                           'tried another arrangement yet (swap or mirror two neighbouring subsystems), try one; otherwise validate the best map and hand it in.')
    return 'go_on', f'NOT FINISHED, do not hand in yet ({len(versions) - 1} of {max_specs} changes used; best so far: {where}). Fix first: {fix}'
