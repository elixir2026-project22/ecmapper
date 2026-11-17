#!/usr/bin/env python3
"""Test suite: python tests/run_tests.py   (offline; the tests of the layout commands need escher_handlayout_kit)"""
import json
import os
import subprocess
import sys
import tempfile
import unittest

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
from escher_rules import probes as P, report, ruleset  # noqa: E402
from escher_rules.mapmodel import MapModel  # noqa: E402

RULES = os.path.join(ROOT, 'rules', 'metabolic_map_rules.yaml')
ENGINE = os.path.join(os.path.dirname(ROOT), 'escher_handlayout_kit')
DRAFT = os.path.join(ENGINE, 'examples', 'central_metabolism_escher.json')
FINAL = os.path.join(ENGINE, 'examples', 'expected_output', 'central_metabolism_escher_final_packed.json')
HAVE_ENGINE = os.path.isfile(DRAFT)


def make_map(nodes, reactions, titles=(), path=None):
    """A small Escher map. nodes: {id: (bigg_id, x, y, primary)}; reactions: {id: (bigg_id, substrate node ids,
    product node ids, reversible)}; the markers lie between the first substrate and the first product, or at the
    point given as a fifth entry. titles: (text, x, y)."""
    N = {n: dict(node_type='metabolite', bigg_id=b, name=b, x=x, y=y, label_x=x + 20, label_y=y, node_is_primary=prim)
         for n, (b, x, y, prim) in nodes.items()}
    R, k = {}, 1000
    for rid, (bigg, subs, prods, rev, *mid) in reactions.items():
        a = nodes[subs[0]] if subs else nodes[prods[0]]
        b = nodes[prods[0]] if prods else (a[0], a[1] + 100, a[2] + 100, True)
        mx, my = mid[0] if mid else ((a[1] + b[1]) / 2, (a[2] + b[2]) / 2)
        m1, mm, m2 = str(k), str(k + 1), str(k + 2)
        k += 3
        N[m1] = dict(node_type='multimarker', x=mx + (a[1] - mx) * 0.2, y=my + (a[2] - my) * 0.2)
        N[mm] = dict(node_type='midmarker', x=mx, y=my)
        N[m2] = dict(node_type='multimarker', x=mx + (b[1] - mx) * 0.2, y=my + (b[2] - my) * 0.2)
        segs = [(m1, mm), (mm, m2)] + [(s, m1) for s in subs] + [(m2, p) for p in prods]
        R[rid] = dict(bigg_id=bigg, name=bigg, reversibility=rev, label_x=mx, label_y=my, gene_reaction_rule='', genes=[],
                      metabolites=[dict(bigg_id=nodes[s][0], coefficient=-1.0) for s in subs]
                      + [dict(bigg_id=nodes[p][0], coefficient=1.0) for p in prods],
                      segments={str(k + i): dict(from_node_id=f, to_node_id=t, b1=None, b2=None) for i, (f, t) in enumerate(segs)})
        k += len(segs)
    T = {str(9000 + i): dict(text=t, x=x, y=y) for i, (t, x, y) in enumerate(titles)}
    data = [dict(map_name='test', map_id='t', map_description='', homepage='', schema=''),
            dict(reactions=R, nodes=N, text_labels=T, canvas=dict(x=-100, y=-100, width=3000, height=3000))]
    if path:
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data, f)
    return data


class Maps(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.addCleanup(self.tmp.cleanup)
        self.rs = ruleset.load([RULES], probes=P.PROBES)
        self.n = 0

    def path(self, name):
        return os.path.join(self.tmp.name, name)

    def model(self, nodes, reactions, titles=(), draft=None):
        self.n += 1
        p = self.path(f'map{self.n}.json')
        make_map(nodes, reactions, titles, p)
        return MapModel(p, self.rs.settings, units_from=draft)

    def run_probe(self, name, model, **params):
        return P.run(name, model, params)


class TestRuleSet(unittest.TestCase):
    def test_bundled_rule_set_is_valid(self):
        rs = ruleset.load(probes=P.PROBES)
        self.assertEqual(rs.problems, [])
        self.assertEqual(len(rs), 66)
        self.assertEqual({r.type for r in rs}, {'motif', 'arrangement', 'encoding'})
        self.assertEqual(sum(1 for r in rs if r.probe), 25)
        self.assertEqual(sorted(r.id for r in rs if r.scope == 'draft'), ['A02', 'M01', 'M02', 'M03'])
        self.assertEqual(rs.get('A05').report_as, 'should')
        self.assertEqual(rs.get('A05').strength, 'must')

    def test_markdown_companion_agrees_with_yaml(self):
        rs = ruleset.load([RULES])
        with open(os.path.join(ROOT, 'rules', 'metabolic_map_rules.md'), encoding='utf-8') as f:
            text = f.read()
        for r in rs:
            self.assertIn(ruleset.rule_row(r), text, msg=r.id)

    def test_render_has_every_rule_and_the_references(self):
        rs = ruleset.load(probes=P.PROBES)
        text = ruleset.to_markdown(rs, bindings=True)
        for r in rs:
            self.assertIn(f'| {r.id} | ', text)
        self.assertIn('`antiparallel_one_axis`', text)
        self.assertIn('## Conflicts between sources', text)
        self.assertIn('32. Purchase H (1997)', text)
        only = ruleset.to_markdown(rs, rs.select(only=['M01']))
        self.assertIn('17. ', only)
        self.assertNotIn('32. Purchase', only)

    def test_select(self):
        rs = ruleset.load()
        self.assertEqual(len(rs.select(strengths=['must'])), 12)
        self.assertEqual([r.id for r in rs.select(groups=['duplication'])], ['E10', 'E11', 'E12', 'E13', 'E14'])
        self.assertEqual([r.id for r in rs.select(only=['M01', 'E10'], skip=['E10'])], ['M01'])
        self.assertEqual(ruleset.expand_ids('A04, E10-E14', {r.id for r in rs}), ['A04', 'E10', 'E11', 'E12', 'E13', 'E14'])

    def _files(self, d, rules, checks=None):
        import yaml
        p = os.path.join(d, 'r.yaml')
        with open(p, 'w', encoding='utf-8') as f:
            yaml.safe_dump(dict(rules=rules, references={1: 'a paper'}), f)
        if checks is not None:
            with open(os.path.join(d, 'r.checks.yaml'), 'w', encoding='utf-8') as f:
                yaml.safe_dump(checks, f)
        return p

    def test_problems_are_all_reported(self):
        good = dict(id='X1', type='motif', when='w', do='d', strength='must', sources=[1])
        with tempfile.TemporaryDirectory() as d:
            p = self._files(d, [good, dict(good), dict(good, id='X2', strength='always', sources=[7]),
                                dict(id='X3', type='motif', when='w', strength='may', colour='red')],
                            dict(checks={'X1': dict(probe='no_such_probe'), 'X2': dict(probe='drawn_once', max_copies=2, report_as='never'),
                                         'X9': dict(probe='drawn_once'), 'X3': dict(scope='elsewhere')}))
            rs = ruleset.load([p], probes=P.PROBES, strict=False)
            text = '\n'.join(rs.problems)
            for part in ('X1: the id is used twice', 'strength `always`', 'source [7]', 'X3: `do` is missing', 'unknown key `colour`',
                         'unknown probe `no_such_probe`', 'no parameter `max_copies`', 'report_as `never`', 'rule X9 is not in the rule set',
                         'scope `elsewhere`'):
                self.assertIn(part, text)
            with self.assertRaises(ruleset.RuleError):
                ruleset.load([p], probes=P.PROBES)
            with self.assertRaises(ruleset.RuleError):
                ruleset.load([os.path.join(d, 'missing.yaml')])

    def test_a_rule_set_without_schema_or_checks_file(self):
        with tempfile.TemporaryDirectory() as d:
            p = self._files(d, [dict(id='X1', type='motif', when='w', do='d', strength='should',
                                     check=dict(probe='currency_not_shared'))])
            rs = ruleset.load([p], probes=P.PROBES)
            self.assertEqual(rs.get('X1').probe, 'currency_not_shared')


class TestProbes(Maps):
    def test_currency_shared_and_metabolite_drawn_twice(self):
        nodes = {'1': ('glc_c', 0, 0, True), '2': ('g6p_c', 0, 400, True), '3': ('f6p_c', 0, 800, True), '4': ('atp_c', 150, 150, False),
                 '5': ('g6p_c', 600, 0, True), '6': ('6pgl_c', 600, 400, True)}
        m = self.model(nodes, {'10': ('HEX1', ['1', '4'], ['2'], False), '11': ('PGI', ['2', '4'], ['3'], True),
                               '12': ('G6PDH', ['5'], ['6'], False)})
        f = self.run_probe('currency_not_shared', m)
        self.assertEqual(len(f.violations), 1)
        self.assertIn('atp_c (node 4)', f.violations[0])
        f = self.run_probe('drawn_once', m)
        self.assertEqual(len(f.violations), 1)
        self.assertIn('g6p_c is drawn 2 times (nodes 2, 5', f.violations[0])
        self.assertEqual(self.run_probe('drawn_once', m, max_nodes=2).violations, [])

    def test_overlap_and_line_through_node(self):
        nodes = {'1': ('a_c', 0, 0, True), '2': ('b_c', 0, 600, True), '3': ('c_c', 5, 300, True), '4': ('d_c', 400, 300, True),
                 '5': ('e_c', 410, 310, True)}
        m = self.model(nodes, {'10': ('R1', ['1'], ['2'], False), '11': ('R2', ['3'], ['4'], False), '12': ('R3', ['4'], ['5'], False)})
        over = self.run_probe('no_node_overlap', m).violations
        self.assertEqual(len(over), 1)
        self.assertIn('d_c (node 4) and e_c (node 5)', over[0])
        through = self.run_probe('edges_avoid_nodes', m).violations
        self.assertTrue(any('R1' in v and 'c_c (node 3)' in v for v in through), through)

    def test_cofactor_pair_and_side_branch(self):
        nodes = {'1': ('glc_c', 0, 0, True), '2': ('g6p_c', 0, 400, True), '3': ('atp_c', 150, 100, False), '4': ('adp_c', -150, 300, False),
                 '5': ('h2o_c', 900, 200, False)}
        m = self.model(nodes, {'10': ('HEX1', ['1', '3', '5'], ['2', '4'], False)})
        f = self.run_probe('cofactor_pairs_same_side', m)
        self.assertEqual((len(f.applies), len(f.violations)), (1, 1))
        self.assertTrue(any('h2o_c (node 5) hangs' in v for v in self.run_probe('side_branches', m).violations))
        nodes['4'] = ('adp_c', 150, 300, False)
        self.assertEqual(self.run_probe('cofactor_pairs_same_side', self.model(nodes, {'10': ('HEX1', ['1', '3'], ['2', '4'], False)})).violations, [])

    def test_antiparallel_reactions(self):
        nodes = {'1': ('f6p_c', 0, 0, True), '2': ('fdp_c', 0, 400, True), '3': ('fdp_c', 500, 400, True)}
        one_axis = self.model(nodes, {'10': ('PFK', ['1'], ['2'], False, (40, 200)), '11': ('FBP', ['2'], ['1'], False, (-40, 200))})
        f = self.run_probe('antiparallel_one_axis', one_axis)
        self.assertEqual((len(f.applies), f.violations), (1, []))
        cycle = self.model(nodes, {'10': ('PFK', ['1'], ['2'], False, (500, 200)), '11': ('FBP', ['2'], ['1'], False, (-40, 200))})
        self.assertIn('two-node cycle', self.run_probe('antiparallel_one_axis', cycle).violations[0])
        copies = self.model(nodes, {'10': ('PFK', ['1'], ['2'], False), '11': ('FBP', ['3'], ['1'], False)})
        self.assertIn('different copies', self.run_probe('antiparallel_one_axis', copies).violations[0])

    def _cycle(self, pts):
        nodes = {str(i + 1): (f'm{i}_c', x, y, True) for i, (x, y) in enumerate(pts)}
        n = len(pts)
        rx = {str(20 + i): (f'R{i}', [str(i + 1)], [str((i + 1) % n + 1)], False) for i in range(n)}
        return self.model(nodes, rx)

    def test_cycle_as_ring(self):
        import math
        ring = self._cycle([(500 * math.cos(2 * math.pi * i / 6), 500 * math.sin(2 * math.pi * i / 6)) for i in range(6)])
        f = self.run_probe('cycles_as_rings', ring)
        self.assertEqual((len(f.applies), f.violations), (1, []))
        flat = self._cycle([(0, 0), (300, 40), (600, 0), (600, 60), (300, 20), (0, 60)])
        self.assertIn('not a ring', self.run_probe('cycles_as_rings', flat).violations[0])
        self.assertEqual(self.run_probe('cycles_as_rings', flat, min_length=7).applies, [])
        chain = self.model({'1': ('a_c', 0, 0, True), '2': ('b_c', 0, 300, True), '3': ('c_c', 0, 600, True)},
                           {'10': ('R1', ['1'], ['2'], False), '11': ('R2', ['2'], ['3'], False)})
        self.assertEqual(self.run_probe('cycles_as_rings', chain).applies, [])

    def test_ring_not_closed(self):
        nodes = {'1': ('a_c', 0, 0, True), '2': ('b_c', 400, 0, True), '3': ('c_c', 200, 350, True), '4': ('a_c', 900, 0, True)}
        m = self.model(nodes, {'10': ('R1', ['1'], ['2'], False), '11': ('R2', ['2'], ['3'], False), '12': ('R3', ['3'], ['4'], False)})
        self.assertIn('not closed, a_c', self.run_probe('cycles_as_rings', m).violations[0])

    def _two_units(self, second_origin):
        ox, oy = second_origin
        nodes = {'1': ('glc_c', 0, 0, True), '2': ('g6p_c', 0, 400, True), '3': ('pyr_c', 0, 800, True),
                 '4': ('pyr_c', ox, oy, True), '5': ('cit_c', ox, oy + 400, True), '6': ('g6p_c', ox + 300, oy, True),
                 '7': ('akg_c', ox + 300, oy + 400, True)}
        rx = {'10': ('R1', ['1'], ['2'], False), '11': ('R2', ['2'], ['3'], False), '12': ('R3', ['4'], ['5'], False),
              '13': ('R4', ['6'], ['7'], False)}
        return self.model(nodes, rx, titles=[('glycolysis', -50, -80), ('citric acid cycle', ox - 50, oy - 80)])

    def test_units_and_their_relations(self):
        m = self._two_units((0, 5000))
        self.assertEqual({u: len(rs) for u, rs in m.units.items()}, {'glycolysis': 2, 'citric acid cycle': 2})
        rel = [dict(unit='glycolysis', **{'is': 'above'}, other='citric acid cycle|TCA')]
        self.assertEqual(self.run_probe('relative_position', m, relations=rel).violations, [])
        far = self.run_probe('linked_units_adjacent', m)
        self.assertEqual(len(far.violations), 1)
        self.assertIn('share 2 metabolites', far.violations[0])
        self.assertEqual(len(self.run_probe('interface_metabolites', m).applies), 2)
        near = self._two_units((900, -3000))
        self.assertEqual(len(self.run_probe('relative_position', near, relations=rel).violations), 1)
        self.assertEqual(self.run_probe('relative_position', near, relations=[dict(unit='lipid', **{'is': 'above'}, other='x')]).applies, [])
        g6p = [dict(unit='citric', **{'is': 'near'}, metabolite='g6p', max_distance=1000)]
        self.assertEqual(len(self.run_probe('relative_position', m, relations=g6p).violations), 1)

    def test_flow_alignment_and_size(self):
        nodes = {'1': ('a_c', 0, 0, True), '2': ('b_c', 0, 400, True), '3': ('c_c', 0, 800, True), '4': ('d_c', 300, 500, True)}
        m = self.model(nodes, {'10': ('R1', ['1'], ['2'], False), '11': ('R2', ['2'], ['3'], False), '12': ('R3', ['3'], ['4'], False)})
        self.assertEqual(self.run_probe('flow_direction', m).violations, [])
        self.assertEqual(len(self.run_probe('flow_direction', m, direction='up').violations), 1)
        self.assertEqual(self.run_probe('axis_aligned', m).violations, [])
        self.assertEqual(len(self.run_probe('axis_aligned', m, min_fraction=0.9).violations), 1)
        self.assertEqual(len(self.run_probe('unit_size', m, max_reactions=2).violations), 1)
        self.assertEqual(len(self.run_probe('unit_size', m, max_map=2, max_reactions=9).violations), 1)
        self.assertEqual(len(self.run_probe('label_legible', m).violations), 1)
        self.assertEqual(self.run_probe('label_legible', m, target_width_mm=2000).violations, [])

    def test_boundary_reactions_and_primary_pair(self):
        nodes = {'1': ('a_c', 1000, 1000, True), '2': ('b_c', 0, 0, True), '3': ('c_c', 2000, 2000, True), '4': ('a_e', 1000, 1200, True),
                 '5': ('coa_c', 300, 300, True)}
        m = self.model(nodes, {'10': ('EX_a_e', ['4'], [], False), '11': ('BIOMASS_x', ['1', '2'], ['3'], False),
                               '12': ('R1', ['2', '5'], ['3'], False), '13': ('TR', ['4'], ['1'], True)})
        f = self.run_probe('boundary_reactions', m)
        self.assertEqual(len(f.applies), 2)
        self.assertTrue(any('BIOMASS_x' in v and 'own copies' in v for v in f.violations))
        self.assertTrue(any('EX_a_e' in v and 'interior' in v for v in f.violations))
        pairs = self.run_probe('one_primary_pair', m).violations
        self.assertTrue(any('R1' in v and '2 primary substrates' in v for v in pairs) and not any('EX_a_e' in v for v in pairs))
        self.assertIn('coa_c (node 5) is primary', self.run_probe('primary_not_carrier', m).violations[0])
        self.assertEqual(len(self.run_probe('transport_reactions', m).applies), 1)
        self.assertIn('2 compartments', self.run_probe('compartments', m).applies[0])

    def test_positions_stable(self):
        nodes = {'1': ('a_c', 0, 0, True), '2': ('b_c', 0, 400, True)}
        rx = {'10': ('R1', ['1'], ['2'], False)}
        old = self.model(nodes, rx)
        new = self.model(dict(nodes, **{'2': ('b_c', 50, 400, True)}), rx)
        self.assertFalse(P.run('positions_stable', new).checked)
        moved = P.run('positions_stable', new, previous=old).violations
        self.assertTrue(any('b_c (node 2) moved 50' in v for v in moved))
        self.assertEqual(P.run('positions_stable', old, previous=old).violations, [])

    def test_every_probe_has_a_test_or_runs(self):
        m = self._two_units((0, 5000))
        for name in P.PROBES:
            self.assertIsInstance(P.run(name, m), P.Finding, msg=name)


class TestReports(Maps):
    def _map(self):
        nodes = {'1': ('glc_c', 0, 0, True), '2': ('g6p_c', 0, 400, True), '3': ('f6p_c', 0, 800, True), '4': ('atp_c', 150, 150, False)}
        p = self.path('m.json')
        make_map(nodes, {'10': ('HEX1', ['1', '4'], ['2'], False), '11': ('PGI', ['2', '4'], ['3'], True)}, [('glycolysis', 0, -80)], p)
        return p

    def test_check_verdict_follows_the_strength_in_the_yaml(self):
        p = self._map()
        _, data, summary = report.check(self.rs, p, out_dir=self.path('out'))
        self.assertEqual(data['state'], 'not_met')
        self.assertEqual(data['broken']['must'], ['E10'])
        self.assertTrue(summary[-1].startswith('RULES NOT MET: 1 must rule is broken (E10)'))
        with open(self.path('out/rules_report.md'), encoding='utf-8') as f:
            text = f.read()
        self.assertIn('atp_c (node 4) is shared by HEX1 (reaction 10), PGI (reaction 11)', text)
        self.assertIn('## To judge by eye', text)
        self.assertIn('## Outside the layout of this map', text)
        # the same rule set with E10 as a should rule: edited as a user would
        with open(RULES, encoding='utf-8') as f:
            rules = f.read()
        edited = rules.replace("never draw it as a shared hub.\n  strength: must", "never draw it as a shared hub.\n  strength: should")
        self.assertNotEqual(edited, rules)
        with open(self.path('mine.yaml'), 'w', encoding='utf-8') as f:
            f.write(edited)
        with open(RULES.replace('.yaml', '.checks.yaml'), encoding='utf-8') as f, open(self.path('mine.checks.yaml'), 'w', encoding='utf-8') as g:
            g.write(f.read())
        rs = ruleset.load([self.path('mine.yaml')], probes=P.PROBES)
        _, data, summary = report.check(rs, p)
        self.assertEqual(data['state'], 'should_open')
        self.assertTrue(summary[-1].startswith('MUST RULES MET'))
        _, data, _ = report.check(self.rs, p, rules=self.rs.select(skip=['E10', 'M23', 'A12']))
        self.assertNotIn('E10', data['broken']['must'])

    def test_brief(self):
        text, data = report.brief(self.rs, self._map(), self.path('b'))
        self.assertTrue(os.path.isfile(self.path('b/rules_brief.md')) and os.path.isfile(self.path('b/rules_brief.json')))
        for part in ('## 1. Rules with places found on this map', '### E10 (must) Currency metabolite', '## 2. Decisions to take once',
                     '**A04**', '**C1**', '## 3. Given by the draft', '## 4. Rules the kit cannot place', '| M13 |'):
            self.assertIn(part, text)
        self.assertEqual(len(data['rules']), 66)

    @unittest.skipUnless(HAVE_ENGINE, 'escher_handlayout_kit is not beside this kit')
    def test_bundled_example(self):
        text, _ = report.brief(self.rs, DRAFT)
        self.assertIn('PFK', text)
        self.assertIn('unit citric acid cycle: acon_C_c - cit_c - oaa_c - mal__L_c - fum_c - succ_c - succoa_c - akg_c - icit_c', text)
        self.assertLess(len(text), 20000)  # fits into one read of a chat workspace
        _, draft, _ = report.check(self.rs, DRAFT)
        _, final, _ = report.check(self.rs, FINAL, draft=DRAFT)
        count = lambda d, rule: len(next(r for r in d['rules'] if r['id'] == rule).get('violations', []))
        self.assertGreater(count(draft, 'E11'), count(final, 'E11'))  # the layout merged shared metabolites
        self.assertGreater(count(draft, 'A17'), 0)
        self.assertEqual(count(final, 'A17'), 0)  # and put linked units next to each other
        self.assertEqual(count(final, 'M06'), 0)
        units = {r.rid: r.unit for r in MapModel(FINAL, self.rs.settings, units_from=DRAFT).reactions.values()}
        self.assertEqual(units, {r.rid: r.unit for r in MapModel(DRAFT, self.rs.settings).reactions.values()})


class TestCommandLine(unittest.TestCase):
    def run_kit(self, *args, cwd=None):
        done = subprocess.run([sys.executable, os.path.join(ROOT, 'run.py'), *args], capture_output=True, text=True, cwd=cwd,
                              encoding='utf-8', env=dict(os.environ, PYTHONIOENCODING='utf-8'))
        return done.returncode, done.stdout + done.stderr

    def test_rule_commands(self):
        code, out = self.run_kit('rules', '--validate')
        self.assertEqual(code, 0, out)
        self.assertIn('rule set ok: 66 rules (25 with a probe)', out)
        code, out = self.run_kit('rules', '--strength', 'must', '--bindings')
        self.assertIn('12 of 66 rules', out)
        self.assertIn('[currency_not_shared]', out)
        code, out = self.run_kit('probes')
        self.assertIn('antiparallel_one_axis  [M07]', out)
        self.assertIn('max_offset=350', out)
        code, out = self.run_kit('render', '--id', 'M07')
        self.assertIn('| M07 | Two antiparallel reactions', out)
        with tempfile.TemporaryDirectory() as d:
            bad = os.path.join(d, 'bad.yaml')
            with open(bad, 'w') as f:
                f.write('rules:\n- {id: X1, type: motif, when: w, do: d, strength: perhaps}\n')
            code, out = self.run_kit('rules', '--validate', '--rules', bad)
            self.assertEqual(code, 1)
            self.assertIn('strength `perhaps`', out)
            code, out = self.run_kit('rules', '--rules', bad)
            self.assertIn('RULE ERROR', out)
            make_map({'1': ('a_c', 0, 0, True), '2': ('b_c', 0, 400, True), '3': ('atp_c', 100, 100, False)},
                     {'10': ('R1', ['1', '3'], ['2'], False), '11': ('R2', ['2', '3'], ['1'], False)}, path=os.path.join(d, 'm.json'))
            code, out = self.run_kit('check', '--input', os.path.join(d, 'm.json'), '--strict', '--out', os.path.join(d, 'o'))
            self.assertEqual(code, 1, out)
            self.assertIn('RULES NOT MET', out)
            code, out = self.run_kit('check', '--input', os.path.join(d, 'm.json'), '--strict', '--skip', 'E10')
            self.assertEqual(code, 0, out)

    @unittest.skipUnless(HAVE_ENGINE, 'escher_handlayout_kit is not beside this kit')
    def test_layout_commands_run_the_engine_and_add_the_rules(self):
        with tempfile.TemporaryDirectory() as d:
            code, out = self.run_kit('inspect', '--input', DRAFT, '--out', 'analysis', cwd=d)
            self.assertEqual(code, 0, out)
            self.assertTrue(os.path.isfile(os.path.join(d, 'analysis', 'analysis.md')))
            self.assertTrue(os.path.isfile(os.path.join(d, 'analysis', 'rules_brief.md')))
            code, out = self.run_kit('start', '--input', DRAFT, '--out', 'spec.json', cwd=d)
            self.assertEqual(code, 0, out)
            code, out = self.run_kit('build', '--input', DRAFT, '--spec', 'spec.json', '--out', 'build_0', '--no-png', cwd=d)
            self.assertEqual(code, 0, out)
            lines = out.strip().splitlines()
            self.assertTrue(any(line.startswith('score') or 'NOT FINISHED' in line or 'TARGETS MET' in line for line in lines), out)
            self.assertTrue(any(line.startswith('rules: ') for line in lines), out)
            self.assertRegex(lines[-1], r'^(RULES NOT MET|MUST RULES MET|ALL CHECKED RULES MET)')
            with open(os.path.join(d, 'build_0', 'rules_report.json'), encoding='utf-8') as f:
                data = json.load(f)
            self.assertEqual(data['draft'], DRAFT)
            self.assertNotIn('M01', data['broken']['must'])  # given by the draft: not counted
            code, out = self.run_kit('build', '--input', DRAFT, '--spec', 'spec.json', '--out', 'build_1', '--no-png', '--no-rules',
                                     '--no-history', cwd=d)
            self.assertNotIn('rules: ', out)
            code, out = self.run_kit('doc', 'format_reference')
            self.assertIn('"merges"', out)
            code, out = self.run_kit('doc', 'guide', '--part', '2')
            self.assertIn('[guide: part 2 of 2]', out)


if __name__ == '__main__':
    unittest.main(verbosity=1)
