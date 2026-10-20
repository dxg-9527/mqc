"""Include biomass rescue nutrients in MQC's nominal 60 mmol C budget.

This configures a diagnostic medium, not a measured cultivation condition or
an experimentally fitted growth rate. Only models with recorded auxiliary
carbon supply are affected. Bounds, not transient solver constraints, carry
the allocation into SBML. Artificial precursor sinks remain explicit.
The budget uses explicit formula carbon atoms; R/X groups are not assigned
invented compositions. Unclassified legacy inputs (including bookkeeping
species) are reported, not silently certified as zero-carbon chemicals.
"""
import math

from mqc.solver_utils import robust_pfba as pfba
from mqc.control.boundary_diagnostics import unresolved_formula


CARBON_BUDGET = 60.0


def input_formula_coverage(model):
    partial, unclassified, active = {}, [], set()
    for reaction in model.reactions:
        if len(reaction.metabolites) != 1:
            continue
        met, c = next(iter(reaction.metabolites.items()))
        if c * (reaction.lower_bound if c < 0 else reaction.upper_bound) <= 0:
            continue
        active.add(reaction.id)
        if unresolved_formula(met):
            if (met.elements or {}).get("C", 0) > 0:
                partial[reaction.id] = met.formula
            else:
                unclassified.append(reaction.id)
    return partial, sorted(unclassified), active


def record_formula_coverage(model, report):
    partial, unclassified, active = input_formula_coverage(model)
    report.update(formula_coverage_version=1, partial_formula_inputs=partial,
                  unclassified_existing_inputs=unclassified,
                  artificial_inputs=sorted(set(report["recorded_sink_reactions"]) & active),
                  carbon_accounting=("explicit_atoms_with_unresolved_input_formulas"
                                     if partial or unclassified else "complete_input_formulas"))


def carbon_inputs(model):
    """Read input direction from stoichiometry, including sink/positive inputs."""
    rows, unknown = [], []
    for reaction in model.reactions:
        if len(reaction.metabolites) != 1:
            continue
        met, coefficient = next(iter(reaction.metabolites.items()))
        capacity = max(0., coefficient * (reaction.lower_bound if coefficient < 0
                                         else reaction.upper_bound))
        if capacity <= 0:
            continue
        elements = met.elements
        if not met.formula or not elements:
            unknown.append(reaction.id)
            continue
        carbon = elements.get("C", 0.)
        if carbon > 0:
            rows.append(dict(reaction_id=reaction.id, metabolite_id=met.id,
                             coefficient=coefficient, carbon=carbon,
                             capacity=capacity, bounds=list(reaction.bounds)))
    return rows, unknown


def validate_supply_budget(model, report, tolerance=1e-6):
    """Check explicit-carbon bounds; do not certify unresolved formula parts."""
    if report.get("status") not in ("allocated", "already_within_budget"):
        return
    rows, unknown = carbon_inputs(model)
    capacity = sum(row["carbon"] * row["capacity"] for row in rows)
    if unknown or capacity > report["carbon_budget"] + tolerance:
        raise RuntimeError("Auxiliary supply budget failed final-bound validation")
    if report.get("formula_coverage_version") == 1:
        partial, unclassified, _ = input_formula_coverage(model)
        if (set(unclassified) & set(report["recorded_auxiliary_reactions"])
                or not set(unclassified) <= set(report["unclassified_existing_inputs"])
                or (report["carbon_accounting"] == "complete_input_formulas"
                    and (partial or unclassified))):
            raise RuntimeError("Auxiliary input formula coverage changed after allocation")


def allocate_auxiliary_supply(model, model_info, model_control_info):
    """Apply MQC's pFBA molecular-fraction allocation to *all* rescue carbon.

    As in get_carbon_proportion, fraction_i = uptake_i / sum(uptake);
    carbon_i = 60 * fraction_i; molecular_bound_i = carbon_i / C_i.
    Unlike the old rescue branch this includes newly added sinks, does not
    require an EX prefix or a prior carbon_source metadata tag, and never
    increases an existing bound. Failed allocation is rolled back and raises;
    an unlimited or zero-growth model cannot silently pass this check.
    """
    auxiliary = set(model_info.get("no_gap_rxn_input", []))
    auxiliary.update(model_info.get("sink_rxn", []))
    rows, unknown = carbon_inputs(model)
    _, unclassified, _ = input_formula_coverage(model)
    triggered = sorted(auxiliary & {row["reaction_id"] for row in rows})
    report = dict(status="not_applicable", carbon_budget=CARBON_BUDGET,
                  budget_basis="explicit_formula_carbon_atoms",
                  trigger_reactions=triggered,
                  recorded_auxiliary_reactions=sorted(auxiliary),
                  recorded_sink_reactions=sorted(set(model_info.get("sink_rxn", []))),
                  unknown_formula_inputs=unknown,
                  experimental_values_used=False, changes=[])
    record_formula_coverage(model, report)
    model_info["auxiliary_supply_budget"] = report
    model_control_info["auxiliary_supply_budget"] = report
    if not triggered and not auxiliary.intersection(set(unknown) | set(unclassified)):
        report["carbon_accounting"] = "not_assessed_no_auxiliary_carbon"
        return report
    report["status"] = "failed"
    unbudgetable = set(unknown) | (auxiliary & set(unclassified))
    if unbudgetable:
        report["reason"] = "cannot_budget_inputs_with_unknown_formula"
        raise RuntimeError(report["reason"] + ": " + ", ".join(sorted(unbudgetable)))
    total_capacity = sum(row["carbon"] * row["capacity"] for row in rows)
    report["carbon_capacity_before"] = total_capacity
    if total_capacity <= CARBON_BUDGET + 1e-6:
        report["status"] = "already_within_budget"
        report["carbon_capacity_after"] = total_capacity
        return report
    original = {row["reaction_id"]: tuple(row["bounds"]) for row in rows}
    try:
        before = model.slim_optimize(error_value=None)
        if not math.isfinite(before) or before <= 1e-7:
            raise RuntimeError("no_positive_growth_before_auxiliary_allocation")
        solution = pfba(model)
        uptake = {row["reaction_id"]: max(0., row["coefficient"] *
                   float(solution.fluxes[row["reaction_id"]])) for row in rows}
        total = sum(uptake.values())
        if not math.isfinite(total) or total <= 0:
            raise RuntimeError("no_carbon_uptake_witness_for_auxiliary_allocation")
        report["growth_before"] = float(before)
        report["actual_carbon_before"] = sum(
            row["carbon"] * uptake[row["reaction_id"]] for row in rows)
        for row in rows:
            rid, coefficient = row["reaction_id"], row["coefficient"]
            reaction = model.reactions.get_by_id(rid)
            cap = min(row["capacity"], CARBON_BUDGET * uptake[rid] / total / row["carbon"])
            lb, ub = reaction.bounds
            bounds = (max(lb, cap / coefficient), ub) if coefficient < 0 else (
                lb, min(ub, cap / coefficient))
            reaction.bounds = bounds
            if bounds != original[rid]:
                report["changes"].append(dict(reaction_id=rid, old_bounds=list(original[rid]),
                    new_bounds=list(bounds), pfba_material_uptake=uptake[rid],
                    carbon_atoms=row["carbon"]))
        after = model.slim_optimize(error_value=None)
        if not math.isfinite(after) or after <= 1e-7:
            raise RuntimeError("no_positive_growth_under_auxiliary_carbon_budget")
        report.update(status="allocated", growth_after=float(after),
                      allocation="original_pfba_molecular_fractions",
                      carbon_capacity_after=sum(row["carbon"] * row["capacity"]
                                                for row in carbon_inputs(model)[0]))
        record_formula_coverage(model, report)
        validate_supply_budget(model, report)
    except Exception as error:
        for rid, bounds in original.items():
            model.reactions.get_by_id(rid).bounds = bounds
        report.update(status="failed", rolled_back=True, reason=str(error))
        raise RuntimeError("Auxiliary carbon supply allocation failed: " + str(error)) from error
    changed = {item["reaction_id"]: item for item in report["changes"]}
    for entry in model_info.get("reactions", []):
        if entry.get("id") in changed:
            entry["bounds"] = tuple(changed[entry["id"]]["new_bounds"])
    # Preserve the original min_set record; explicitly record the *final*
    # supply, including rescue inputs absent from the old metadata table.
    final = carbon_inputs(model)[0]
    model_info["carbon_source_boundary_dict"] = {
        row["reaction_id"]: tuple(row["bounds"]) for row in final}
    model_info["carbon_source_boundary"] = [
        "{} : {}".format(row["reaction_id"], tuple(row["bounds"])) for row in final]
    model_control_info.setdefault("boundary_information", {})[
        "carbon_source_boundary"] = ",".join(model_info["carbon_source_boundary"])
    return report
