"""Budget rescue inputs without changing pathways, objectives or normal media."""
import tempfile
import unittest
from pathlib import Path

import cobra
from mqc.control.auxiliary_supply import (
    allocate_auxiliary_supply, carbon_inputs, validate_supply_budget,
)


def toy(positive=False):
    model = cobra.Model("budget")
    a = cobra.Metabolite("sugar_c", formula="C6H12O6", compartment="c")
    b = cobra.Metabolite("precursor_c", formula="C3H6O3", compartment="c")
    exchange = cobra.Reaction("main", lower_bound=-10, upper_bound=1000)
    exchange.add_metabolites({a: -1})
    supply = cobra.Reaction("rescue", lower_bound=-1000 if not positive else 0,
                            upper_bound=0 if not positive else 500)
    supply.add_metabolites({b: -1 if not positive else 2})
    route = cobra.Reaction("route")
    route.add_metabolites({b: -2, a: 1})
    biomass = cobra.Reaction("bio")
    biomass.add_metabolites({a: -1})
    model.add_reactions([exchange, supply, route, biomass])
    model.objective = biomass
    return model


class AuxiliarySupplyTests(unittest.TestCase):
    def test_budget_and_sbml_and_idempotence(self):
        for positive in (False, True):
            with self.subTest(positive=positive):
                model = toy(positive)
                info = dict(no_gap_rxn_input=["rescue"], reactions=[dict(id="rescue")])
                report = allocate_auxiliary_supply(model, info, {})
                self.assertEqual(report["status"], "allocated")
                self.assertGreater(model.slim_optimize(), 0)
                self.assertLessEqual(sum(r["carbon"] * r["capacity"]
                                         for r in carbon_inputs(model)[0]), 60 + 1e-6)
                self.assertEqual(model.reactions.main.upper_bound, 1000)
                self.assertEqual(model.reactions.route.bounds, (0, 1000))
                with tempfile.TemporaryDirectory() as tmp:
                    filename = str(Path(tmp) / "m.xml")
                    cobra.io.write_sbml_model(model, filename)
                    reloaded = cobra.io.read_sbml_model(filename)
                    validate_supply_budget(reloaded, report)
                    self.assertAlmostEqual(model.slim_optimize(), reloaded.slim_optimize())
                bounds = {r.id: r.bounds for r in model.reactions}
                again = allocate_auxiliary_supply(model, info, {})
                self.assertEqual(again["status"], "already_within_budget")
                self.assertEqual(bounds, {r.id: r.bounds for r in model.reactions})

    def test_sink_without_exchange_prefix_is_budgeted(self):
        model = toy()
        report = allocate_auxiliary_supply(model, dict(sink_rxn=["rescue"]), {})
        self.assertEqual(report["status"], "allocated")

    def test_normal_medium_unchanged(self):
        model = toy()
        bounds = {r.id: r.bounds for r in model.reactions}
        self.assertEqual(allocate_auxiliary_supply(model, {}, {})["status"], "not_applicable")
        self.assertEqual(bounds, {r.id: r.bounds for r in model.reactions})

    def test_export_only_sink_does_not_trigger(self):
        model = toy()
        model.reactions.rescue.bounds = (0, 1000)
        report = allocate_auxiliary_supply(model, dict(sink_rxn=["rescue"]), {})
        self.assertEqual(report["status"], "not_applicable")
        self.assertEqual(report["artificial_inputs"], [])
        self.assertEqual(report["recorded_sink_reactions"], ["rescue"])

    def test_unknown_auxiliary_formula_is_not_free_carbon(self):
        model = toy()
        model.metabolites.precursor_c.formula = None
        info = dict(no_gap_rxn_input=["rescue"])
        with self.assertRaisesRegex(RuntimeError, "unknown_formula"):
            allocate_auxiliary_supply(model, info, {})
        self.assertEqual(info["auxiliary_supply_budget"]["status"], "failed")
        self.assertEqual(model.reactions.rescue.bounds, (-1000, 0))

    def test_failed_allocation_restores_bounds(self):
        model = toy()
        model.reactions.bio.lower_bound = 11  # more C than the whole budget
        before = {r.id: r.bounds for r in model.reactions}
        info = dict(no_gap_rxn_input=["rescue"])
        with self.assertRaises(RuntimeError):
            allocate_auxiliary_supply(model, info, {})
        self.assertTrue(info["auxiliary_supply_budget"]["rolled_back"])
        self.assertEqual(before, {r.id: r.bounds for r in model.reactions})

    def test_guard_detects_reopening_and_never_credits_exports(self):
        model = toy()
        report = allocate_auxiliary_supply(model, dict(no_gap_rxn_input=["rescue"]), {})
        model.reactions.rescue.lower_bound = -1000
        with self.assertRaisesRegex(RuntimeError, "final-bound"):
            validate_supply_budget(model, report)

    def test_no_change_when_already_budgeted(self):
        model = toy()
        model.reactions.main.lower_bound = -5
        model.reactions.rescue.lower_bound = -10
        before = {r.id: r.bounds for r in model.reactions}
        report = allocate_auxiliary_supply(model, dict(no_gap_rxn_input=["rescue"]), {})
        self.assertEqual(report["status"], "already_within_budget")
        self.assertEqual(before, {r.id: r.bounds for r in model.reactions})

    def test_partial_formula_is_explicit_not_certified_complete(self):
        model = toy()
        model.metabolites.precursor_c.formula = "C3H6O3R"
        report = allocate_auxiliary_supply(model, dict(no_gap_rxn_input=["rescue"]), {})
        self.assertEqual(report["status"], "allocated")
        self.assertEqual(report["partial_formula_inputs"], {"rescue": "C3H6O3R"})
        self.assertNotEqual(report["carbon_accounting"], "complete_input_formulas")
        validate_supply_budget(model, report)
        report["carbon_accounting"] = "complete_input_formulas"
        with self.assertRaisesRegex(RuntimeError, "formula coverage"):
            validate_supply_budget(model, report)

    def test_unclassified_new_auxiliary_is_not_assumed_carbon_free(self):
        for formula in ("R", "X"):
            model = toy()
            model.metabolites.precursor_c.formula = formula
            with self.assertRaisesRegex(RuntimeError, "unknown_formula"):
                allocate_auxiliary_supply(model, dict(no_gap_rxn_input=["rescue"]), {})

    def test_legacy_symbolic_input_is_flagged_not_deleted(self):
        model = toy()
        marker = cobra.Metabolite("marker_c", formula="X", compartment="c")
        bookkeeping = cobra.Reaction("bookkeeping", lower_bound=0, upper_bound=1000)
        bookkeeping.add_metabolites({marker: 1})
        model.add_reactions([bookkeeping])
        report = allocate_auxiliary_supply(model, dict(no_gap_rxn_input=["rescue"]), {})
        self.assertEqual(report["status"], "allocated")
        self.assertIn("bookkeeping", report["unclassified_existing_inputs"])
        self.assertEqual(bookkeeping.bounds, (0, 1000))
        self.assertNotEqual(report["carbon_accounting"], "complete_input_formulas")
        validate_supply_budget(model, report)


if __name__ == "__main__":
    unittest.main()
