"""Regressions for failures observed in the 2026-10-19 cluster pilot."""
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO/'cluster'))
import cobra
from run_model import load_exported_model
from mqc.utils import keep_decimal_places
import mqc.control.yield_control as yield_module


class ClusterCompatibilityTests(unittest.TestCase):
    def test_sbml_export_with_hidden_xml_filename(self):
        model = cobra.Model('export_regression')
        metabolite = cobra.Metabolite('a_c', compartment='c')
        source = cobra.Reaction('SOURCE')
        source.add_metabolites({metabolite: 1})
        demand = cobra.Reaction('DEMAND')
        demand.add_metabolites({metabolite: -1})
        model.add_reactions([source, demand])
        model.objective = demand
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder)/'.xml'
            cobra.io.write_sbml_model(model, str(path))
            self.assertEqual(path.suffix, '')
            reopened = load_exported_model(path)
        self.assertEqual(set(reopened.reactions.list_attr('id')), {'SOURCE', 'DEMAND'})
        self.assertEqual(reopened.reactions.DEMAND.objective_coefficient, 1.0)
        self.assertEqual(reopened.metabolites.a_c.compartment, 'c')

    def test_positive_exponents_do_not_create_phantom_metabolites(self):
        for bounds in [(0, 1000), (-1000, 0), (-1000, 1000)]:
            with self.subTest(bounds=bounds):
                model = cobra.Model('coefficient_regression')
                metabolites = [cobra.Metabolite(name, compartment='c')
                               for name in ['a_c', 'b_c', 'd_c', 'product_c']]
                reaction = cobra.Reaction('BIOMASS')
                reaction.add_metabolites(dict(zip(metabolites, [-1.0, -41.29123, -4.1267e-5, 453.512])))
                reaction.bounds = bounds
                model.add_reactions([reaction])
                expected = {met.id: float(format(coefficient, '.3e'))
                            for met, coefficient in reaction.metabolites.items()}
                _, equation = keep_decimal_places(reaction)
                reaction.reaction = equation
                self.assertEqual(set(model.metabolites.list_attr('id')), set(expected))
                self.assertEqual({met.id: value for met, value in reaction.metabolites.items()}, expected)
                self.assertEqual(reaction.bounds, bounds)
                self.assertTrue(all(m.compartment == 'c' for m in model.metabolites))

    def test_empty_revision_branch_generates_report_without_model_changes(self):
        model = cobra.Model('report_regression')
        before = cobra.io.to_json(model)
        info = dict(substance_approach=[], total_carbon_content=0,
                    yield_exceeds_the_theoretical_maximum_number=0, yield_net_matter=[],
                    yield_generation=[], yield_revision=[], control_analysis=[])
        control = {'check_metabolite_yield': {}}
        checker = object.__new__(yield_module.Yields)
        checker.find_yield_generation_multiprocessing = lambda i, m, c: c['check_metabolite_yield'].update(score=0)
        checker.get_associated_substances = lambda *args: None
        with mock.patch.multiple(yield_module, set_c_source_supply=lambda *args: None,
                                 set_auto_source_supply=lambda *args: None, is_autotrophic=lambda *args: False):
            checker.yield_control(info, model, model.copy(), control, None, {})
        self.assertIn('[]', control['check_metabolite_yield']['Details'])
        self.assertEqual(cobra.io.to_json(model), before)


if __name__ == '__main__':
    unittest.main()
