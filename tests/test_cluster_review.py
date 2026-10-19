"""Regression coverage for the 213-model cluster review (no phenotype fitting)."""
import unittest
from unittest import mock
from types import SimpleNamespace

import cobra
import pandas as pd

from mqc.control.boundary_diagnostics import (
    unresolved_boundary_source, unresolved_formula, without_unresolved_inputs,
)
from mqc.control import postgap_control as guard
from mqc.utils import max_reduced_degree_rate


class ReviewTests(unittest.TestCase):
    def test_unknown_material_outlet_is_not_a_precursor_source(self):
        r = cobra.Reaction('not_a_named_biomass_exchange', lower_bound=-1000)
        m = cobra.Metabolite('placeholder', formula='X', compartment='c')
        r.add_metabolites({m: -2})
        self.assertFalse(unresolved_boundary_source(r, 3))
        self.assertTrue(unresolved_boundary_source(r, -3))
        r *= -1
        self.assertFalse(unresolved_boundary_source(r, -3))
        self.assertTrue(unresolved_boundary_source(r, 3))
        m.formula = 'C6H12O6'
        self.assertFalse(unresolved_boundary_source(r, 3))

    def test_rhenium_is_not_an_R_group(self):
        m = cobra.Metabolite('m', formula='Re')
        self.assertFalse(unresolved_formula(m))
        for formula in (None, '', 'R', 'C5H9R', 'X', 'nan'):
            m.formula = formula
            self.assertTrue(unresolved_formula(m))

    def test_unknown_inputs_restore_on_exception_without_undoing_internal_repair(self):
        model = cobra.Model('temporary_assay')
        m = cobra.Metabolite('unknown', formula=None, compartment='c')
        n = cobra.Metabolite('known', formula='CH2O', compartment='c')
        for rid, stoich in [('sink', {m:-1}), ('reversed_sink', {m:2}), ('internal', {m:-1,n:1})]:
            r = cobra.Reaction(rid, lower_bound=-1000)
            r.add_metabolites(stoich); model.add_reactions([r])
        with self.assertRaisesRegex(RuntimeError, 'intentional'):
            with without_unresolved_inputs(model) as changes:
                self.assertEqual(len(changes), 2)
                self.assertEqual(model.reactions.sink.bounds, (0,1000))
                self.assertEqual(model.reactions.reversed_sink.bounds, (-1000,0))
                model.reactions.internal.lower_bound = 0
                raise RuntimeError('intentional')
        self.assertEqual(model.reactions.sink.bounds, (-1000,1000))
        self.assertEqual(model.reactions.reversed_sink.bounds, (-1000,1000))
        self.assertEqual(model.reactions.internal.bounds, (0,1000))

    def test_N_only_output_and_signed_scaled_boundary_accounting(self):
        model = cobra.Model('nitrous_oxide_oxidation')
        mets = {k: cobra.Metabolite(k+'_c', name=k, formula=f, charge=0, compartment='c')
                for k,f in [('glycerol','C3H8O3'),('n2o','N2O'),('n2','N2'),('glyald','C3H6O3')]}
        flux = {}
        for key, coefficient, value in [('glycerol',-1,-20), ('n2o',2,10), ('n2',-3,20/3), ('glyald',-1,20)]:
            r = cobra.Reaction('EX_'+key, lower_bound=-1000)
            r.add_metabolites({mets[key]:coefficient}); model.add_reactions([r]); flux[r.id] = value
        with mock.patch.object(model, 'slim_optimize', return_value=20), \
                mock.patch('mqc.utils.pfba', return_value=SimpleNamespace(fluxes=pd.Series(flux))):
            rate, notes = max_reduced_degree_rate(model, mets['glyald'].id, {})
        self.assertAlmostEqual(rate, 20)
        self.assertEqual(notes, [])

    def test_oxygen_secretion_does_not_inflate_yield_ceiling(self):
        model = cobra.Model('oxygen_not_a_yield_credit')
        flux = {}
        for name, formula, value in [('glc','C6H12O6',-10), ('etoh','C2H6O',20.5), ('o2','O2',1.5)]:
            met = cobra.Metabolite(name+'_c', name='O2' if name=='o2' else name,
                                  formula=formula, charge=0, compartment='c')
            r = cobra.Reaction('EX_'+name, lower_bound=-1000)
            r.add_metabolites({met:-1}); model.add_reactions([r]); flux[r.id] = value
        with mock.patch.object(model, 'slim_optimize', return_value=20.5), \
                mock.patch('mqc.utils.pfba', return_value=SimpleNamespace(fluxes=pd.Series(flux))):
            rate, _ = max_reduced_degree_rate(model, 'etoh_c', {})
        self.assertAlmostEqual(rate, 20)

    def guard_fixture(self):
        model = cobra.Model('parallel_cycles')
        a,b = [cobra.Metabolite(x, compartment='c') for x in ('a','b')]
        for name in ('biomass','new_gap','old1','old2'):
            r = cobra.Reaction(name, lower_bound=-1000 if name.startswith('old') else 0)
            r.add_metabolites({a:-1,b:1}); model.add_reactions([r])
        target = dict(target='ATP@c', carrier='ATP', group='energy')
        def diagnose(m, t, tolerance=1e-7, witness=False):
            opened = [r for r in m.reactions if r.id.startswith('old') and r.lower_bound < 0]
            result = dict(target=t['target'], group='energy', value=len(opened),
                          passed=not opened, solver_status='optimal')
            if witness:
                result['witness'] = [dict(reaction_id='new_gap',flux=1)] + [dict(reaction_id=r.id,flux=-1) for r in opened]
            return result
        def growth(m, biomass):
            return dict(value=1 if m.reactions.new_gap.upper_bound else 0, solver_status='optimal')
        return model,target,diagnose,growth

    def test_sequential_existing_cuts_commit_only_when_all_pass(self):
        model,target,diagnose,growth = self.guard_fixture()
        for budget, expected in [(1,'failed'), (2,'repaired')]:
            with mock.patch.object(guard,'discover_targets',return_value=([target],[])), \
                    mock.patch.object(guard,'diagnose',side_effect=diagnose), \
                    mock.patch.object(guard,'_growth',side_effect=growth):
                result = guard.validate_and_repair(model,'biomass',['new_gap'], {'biomass','old1','old2'},
                                                   max_existing_edits=budget)
            self.assertEqual(result['status'], expected)
            if budget == 1:
                self.assertTrue(result['trial_changes_rolled_back'])
                self.assertEqual(model.reactions.old1.lower_bound, -1000)
                self.assertEqual(model.reactions.old2.lower_bound, -1000)
                self.assertEqual(result['repairs'], [])
            else:
                self.assertEqual(len(result['repairs']), 2)
                self.assertEqual(result['energy_status'], 'passed')
                self.assertEqual(result['growth_status'], 'positive')

    def test_energy_pass_with_zero_growth_is_not_completion(self):
        model,target,_,_ = self.guard_fixture()
        with mock.patch.object(guard,'discover_targets',return_value=([target],[])), \
                mock.patch.object(guard,'diagnose',return_value=dict(target='ATP@c',passed=True,solver_status='optimal',value=0)), \
                mock.patch.object(guard,'_growth',return_value=dict(value=0,solver_status='optimal')):
            result = guard.validate_and_repair(model,'biomass',[],set())
        self.assertEqual(result['status'],'failed')
        self.assertEqual(result['reason'],'zero_growth_after_correction')
        self.assertEqual(result['energy_status'],'passed')
        self.assertEqual(result['growth_status'],'nonpositive_or_unsolved')

    def plateau_fixture(self):
        model, target, diagnose, growth = self.guard_fixture()
        def parallel(m, t, tolerance=1e-7, witness=False):
            result = diagnose(m, t, tolerance, witness)
            result['value'] = 1 if result['value'] else 0
            return result
        return model, target, parallel, growth

    def test_plateau_combo_commits_only_complete_sequence(self):
        model, target, diagnose, growth = self.plateau_fixture()
        with mock.patch.object(guard, 'discover_targets', return_value=([target], [])), \
             mock.patch.object(guard, 'diagnose', side_effect=diagnose), \
             mock.patch.object(guard, '_growth', side_effect=growth):
            result = guard.validate_and_repair(model, 'biomass', ['new_gap'],
                {'biomass','old1','old2'}, max_existing_edits=2)
        self.assertEqual(result['status'], 'repaired')
        self.assertEqual(len(result['repairs']), 2)
        self.assertEqual(model.reactions.new_gap.bounds, (0, 1000))
        self.assertEqual(model.reactions.old1.lower_bound, 0)
        self.assertEqual(model.reactions.old2.lower_bound, 0)
        search = result['attempts'][0]['plateau_search']
        self.assertEqual(search['status'], 'validated_completion')
        self.assertTrue(all(t['completion']['max_plateau_trials'] == 0 for t in search['trials']))

    def test_plateau_disabled_or_budget_insufficient_restores_all_bounds(self):
        for kwargs in [dict(max_plateau_trials=0), dict(max_existing_edits=1), dict(max_edits=1)]:
            model, target, diagnose, growth = self.plateau_fixture()
            before = {r.id:r.bounds for r in model.reactions}
            with mock.patch.object(guard, 'discover_targets', return_value=([target], [])), \
                 mock.patch.object(guard, 'diagnose', side_effect=diagnose), \
                 mock.patch.object(guard, '_growth', side_effect=growth):
                result = guard.validate_and_repair(model, 'biomass', ['new_gap'],
                    {'biomass','old1','old2'}, **kwargs)
            self.assertEqual(result['status'], 'failed')
            self.assertEqual(result['repairs'], [])
            self.assertTrue(result['trial_changes_rolled_back'])
            self.assertEqual(before, {r.id:r.bounds for r in model.reactions})

    def test_plateau_trial_limit_is_respected(self):
        model, target, diagnose, growth = self.plateau_fixture()
        with mock.patch.object(guard, 'discover_targets', return_value=([target], [])), \
             mock.patch.object(guard, 'diagnose', side_effect=diagnose), \
             mock.patch.object(guard, '_growth', side_effect=growth):
            result = guard.validate_and_repair(model, 'biomass', ['new_gap'],
                {'biomass','old1','old2'}, max_plateau_trials=1)
        search = result['attempts'][0]['plateau_search']
        self.assertEqual(search['tested'], 1)
        self.assertEqual(search['untested_candidates'], 1)
        self.assertEqual(result['status'], 'repaired')

    def test_plateau_cannot_reactivate_a_passed_carrier(self):
        model, target, diagnose, growth = self.plateau_fixture()
        redox = dict(target='NADH@c', carrier='NADH', group='redox')
        before = {r.id:r.bounds for r in model.reactions}
        def checks(m, t, tolerance=1e-7, witness=False):
            if t['target'] == 'ATP@c':
                return diagnose(m,t,tolerance,witness)
            activated = any(m.reactions.get_by_id(r).lower_bound == 0 for r in ('old1','old2'))
            return dict(target=t['target'], group='redox', value=int(activated),
                        passed=not activated, solver_status='optimal')
        with mock.patch.object(guard, 'discover_targets', return_value=([target,redox], [])), \
             mock.patch.object(guard, 'diagnose', side_effect=checks), \
             mock.patch.object(guard, '_growth', side_effect=growth):
            result = guard.validate_and_repair(model, 'biomass', ['new_gap'], {'biomass','old1','old2'})
        self.assertEqual(result['status'], 'failed')
        self.assertEqual(result['attempts'][0]['plateau_search']['tested'], 0)
        self.assertEqual(before, {r.id:r.bounds for r in model.reactions})

    def test_plateau_exception_rolls_back_entire_trial(self):
        model, target, diagnose, growth = self.plateau_fixture()
        before = {r.id:r.bounds for r in model.reactions}
        calls = {'count':0}
        def discover(m):
            calls['count'] += 1
            if calls['count'] > 1:
                raise RuntimeError('intentional branch failure')
            return [target], []
        with mock.patch.object(guard, 'discover_targets', side_effect=discover), \
             mock.patch.object(guard, 'diagnose', side_effect=diagnose), \
             mock.patch.object(guard, '_growth', side_effect=growth):
            with self.assertRaisesRegex(RuntimeError, 'intentional branch failure'):
                guard.validate_and_repair(model, 'biomass', ['new_gap'], {'biomass','old1','old2'})
        self.assertEqual(before, {r.id:r.bounds for r in model.reactions})

    def test_missing_biomass_is_not_invented(self):
        model = cobra.Model('no_biomass')
        result = guard.validate_and_repair(model, '', [], set())
        self.assertEqual(result['reason'], 'missing_biomass_objective')
        self.assertEqual(len(model.reactions),0)


if __name__ == '__main__':
    unittest.main()
