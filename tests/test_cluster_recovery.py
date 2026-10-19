"""Regression tests for cluster failures, independent of strain phenotypes."""
import json
import math
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT), str(ROOT / 'cluster')]
import cobra
from cobra.exceptions import Infeasible, OptimizationError
from mqc import solver_utils as numeric
from mqc.control.formula_charge_pairing import zero_charge_placeholders, complete_reference_pair
import common
import run_model
from mqc.control.nitric_oxide_direction import nitric_oxide_direction, nitric_oxide_bounds


def toy():
    model = cobra.Model('numerical_test')
    met = cobra.Metabolite('a_c', formula='CH2O', compartment='c')
    source = cobra.Reaction('SOURCE', upper_bound=10)
    demand = cobra.Reaction('DEMAND')
    source.add_metabolites({met: 1})
    demand.add_metabolites({met: -1})
    model.add_reactions([source, demand])
    model.objective = demand
    model.add_cons_vars(model.problem.Constraint(source.flux_expression, ub=3, name='custom_capacity'))
    return model


class RecoveryTests(unittest.TestCase):
    def test_unrounded_yield(self):
        self.assertFalse(numeric.exceeds_yield_limit(9.375, 9.374998571428572))
        self.assertTrue(numeric.exceeds_yield_limit(9.375, 9.37))
        self.assertTrue(numeric.exceeds_yield_limit(1e-4, 0))
        self.assertFalse(numeric.exceeds_yield_limit(1e-8, 0))
        for value in [math.nan, math.inf]:
            with self.assertRaises(OptimizationError):
                numeric.exceeds_yield_limit(value, 1)

    def test_cold_copy_preserves_custom_constraint(self):
        model = toy()
        copied = numeric.cold_copy(model)
        self.assertEqual(set(model.constraints.keys()), set(copied.constraints.keys()))
        self.assertAlmostEqual(copied.slim_optimize(error_value=None), 3)
        self.assertLessEqual(numeric.validate_primals(copied), 1e-6)
        self.assertEqual(model.reactions.SOURCE.upper_bound, 10)

    def test_cold_retry_only_after_known_failure(self):
        model = toy()
        original = cobra.Model.slim_optimize
        calls = []
        def fail_first(instance, *args, **kwargs):
            calls.append(1)
            if len(calls) == 1:
                raise Infeasible('simulated warm state')
            return original(instance, *args, **kwargs)
        with mock.patch.object(cobra.Model, 'slim_optimize', fail_first):
            cold, value, info = numeric.solve_with_cold_retry(model)
        self.assertAlmostEqual(value, 3)
        self.assertEqual(info['method'], 'complete_lp_cold_retry')
        self.assertIsNot(cold, model)
        with mock.patch.object(model, 'slim_optimize', side_effect=ValueError('programming error')):
            with self.assertRaises(ValueError):
                numeric.solve_with_cold_retry(model)

    def test_truly_infeasible_not_marked_as_pass(self):
        model = toy()
        model.reactions.DEMAND.lower_bound = 4
        with self.assertRaises(OptimizationError):
            numeric.solve_with_cold_retry(model)

    def test_pfba_retry_preserves_caller_and_custom_constraint(self):
        model = toy()
        before = cobra.io.to_json(model)
        original = numeric.cobra_pfba
        calls = []
        def fail_first(*args, **kwargs):
            calls.append(1)
            if len(calls) == 1:
                raise Infeasible('simulated strict pFBA failure')
            return original(*args, **kwargs)
        with mock.patch.object(numeric, 'cobra_pfba', side_effect=fail_first):
            solution = numeric.robust_pfba(model)
        self.assertAlmostEqual(solution.fluxes['DEMAND'], 3)
        self.assertEqual(cobra.io.to_json(model), before)
        self.assertIn('custom_capacity', model.constraints)

    def test_pfba_relaxation_is_bounded(self):
        model = toy()
        original = numeric._checked_pfba
        def fail_exact(work, fraction):
            if fraction == 1.0:
                raise Infeasible('simulated exact optimum instability')
            return original(work, fraction)
        with mock.patch.object(numeric, 'cobra_pfba', side_effect=Infeasible('warm')), \
                mock.patch.object(numeric, '_checked_pfba', side_effect=fail_exact):
            solution = numeric.robust_pfba(model)
        self.assertLessEqual(abs(solution.fluxes['DEMAND'] - 3), 3.001e-6)

    def test_git_permission_failure_uses_content_lock(self):
        with mock.patch.object(common.subprocess, 'run', side_effect=PermissionError('git')):
            self.assertEqual(common.git_commit(), 'unavailable')

    def test_startup_failure_is_recorded(self):
        with tempfile.TemporaryDirectory() as folder:
            argv = ['run_model.py', '--input', str(Path(folder)/'missing.xml'), '--work', folder]
            with mock.patch.object(sys, 'argv', argv), mock.patch.object(run_model, 'verify_release', side_effect=ValueError('bad lock')):
                self.assertEqual(run_model.main(), 1)
            result = json.loads((Path(folder)/'outcome.json').read_text())
            self.assertEqual(result['status'], 'failed')
            self.assertIn('bad lock', result['error'])
            self.assertEqual(result['stages'], [])

    def test_explicit_neutral_not_generally_missing(self):
        model = toy()
        self.assertEqual(zero_charge_placeholders(model), set())
        for name, formula in [('h_c','H'),('na1_c','Na'),('cl_c','Cl')]:
            model.add_metabolites([cobra.Metabolite(name, formula=formula, charge=0, compartment='c')])
        model.metabolites.a_c.charge = 0
        self.assertEqual(zero_charge_placeholders(model), set(model.metabolites.list_attr('id')))
        model.metabolites.h_c.charge = 1
        self.assertEqual(zero_charge_placeholders(model), set())

    def test_placeholder_completion_requires_exact_formula(self):
        met = cobra.Metabolite('a', formula='C2H4O2', charge=0)
        met.notes['MQC_ZERO_CHARGE_PLACEHOLDER'] = 'recognized'
        self.assertFalse(complete_reference_pair(met, 'C2H3O2', -1))
        self.assertEqual(met.charge, 0)
        self.assertTrue(complete_reference_pair(met, 'C2H4O2', 0))

    def test_nitric_oxide_rule_is_stoichiometric_not_reaction_id(self):
        model = cobra.Model('NO_test')
        specs = [('no_e', 'NO', -2), ('n2o_e', 'N2O', 1), ('h2o_e', 'H2O', 1),
                 ('q8_c', 'C49H74O4', 1), ('q8h2_c', 'C49H76O4', -1)]
        reaction = cobra.Reaction('arbitrary_identifier', lower_bound=-1000)
        reaction.annotation['ec-code'] = '1.7.5.2'
        reaction.add_metabolites({cobra.Metabolite(mid, formula=formula, charge=0,
            compartment=mid[-1]): coefficient for mid, formula, coefficient in specs})
        model.add_reactions([reaction])
        self.assertEqual(nitric_oxide_direction(reaction), 1)
        self.assertEqual(nitric_oxide_bounds(reaction), (0, 1000))
        reaction *= -1
        self.assertEqual(nitric_oxide_direction(reaction), -1)
        self.assertEqual(nitric_oxide_bounds(reaction), (-1000, 0))
        reaction.metabolites[next(m for m in reaction.metabolites if m.id=='n2o_e')]
        reaction.add_metabolites({model.metabolites.no_e: 1})
        self.assertIsNone(nitric_oxide_direction(reaction))

    def test_postgap_existing_fallback_is_transactional(self):
        from mqc.control import postgap_control as guard
        model = toy()
        extra = cobra.Reaction('NEW_GAP')
        extra.add_metabolites({model.metabolites.a_c:-1,
            cobra.Metabolite('b_c', compartment='c'):1})
        old = cobra.Reaction('OLD_TRANSPORT', lower_bound=-1000)
        old.add_metabolites({model.metabolites.a_c:-1,
            cobra.Metabolite('a_e', compartment='e'):1})
        model.add_reactions([extra,old])
        target = dict(target='ATP@c', carrier='ATP', group='energy')
        def diagnosis(m, t, tolerance=1e-7, witness=False):
            passed = m.reactions.OLD_TRANSPORT.lower_bound == 0
            result = dict(target='ATP@c', group='energy', value=0 if passed else 1,
                          passed=passed, solver_status='optimal')
            if witness:
                result['witness'] = [dict(reaction_id='NEW_GAP', flux=1),
                                     dict(reaction_id='OLD_TRANSPORT', flux=-1)]
            return result
        def growth(m, biomass):
            return dict(value=0 if m.reactions.NEW_GAP.upper_bound==0 else 1,
                        solver_status='optimal')
        with mock.patch.object(guard, 'discover_targets', return_value=([target],[])), \
                mock.patch.object(guard, 'diagnose', side_effect=diagnosis), \
                mock.patch.object(guard, '_growth', side_effect=growth):
            report = guard.validate_and_repair(model, 'DEMAND', ['NEW_GAP'], {'SOURCE','DEMAND','OLD_TRANSPORT'})
            self.assertEqual(report['status'], 'repaired')
            self.assertEqual(model.reactions.OLD_TRANSPORT.bounds,(0,1000))
            self.assertEqual(model.reactions.NEW_GAP.bounds,(0,1000))
            model.reactions.OLD_TRANSPORT.lower_bound = -1000
            # An existing-only witness is NOT permission to edit native models.
            report = guard.validate_and_repair(model, 'DEMAND', [], set(model.reactions.list_attr('id')))
            self.assertEqual(report['status'], 'failed')
            self.assertEqual(model.reactions.OLD_TRANSPORT.bounds,(-1000,1000))


if __name__ == '__main__':
    unittest.main()
