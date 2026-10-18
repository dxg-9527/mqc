"""Conservative, experiment-independent validation after biomass gap filling.

Checks source-free dissipation of recognized cytosolic energy/redox carriers.
Only newly gap-filled reaction directions may be restricted. Repairs must
preserve growth viability in the unchanged, stored QC medium; an optimum
inflated by a witnessed source-free cycle is not a conservation target. This is
a bounded greedy safeguard, not a thermodynamic proof or a minimum-edit solver.
No experimental growth values or strain-specific reaction lists are used.
"""
import json
import math
from pathlib import Path

from cobra import Reaction
from cobra.flux_analysis import pfba
from mqc.control.growth_environment import validate_growth_environment


class PostGapValidationError(RuntimeError):
    """The final model must not be reported as successfully validated."""


# Exact identifiers take precedence over conservative, exact name aliases.
# Deliberately do not use legacy NADH aliases that also contain NADP(+).
ALIASES = {
    "atp": {"atp", "cpd00002"}, "adp": {"adp", "cpd00008"},
    "h2o": {"h2o", "cpd00001"}, "h": {"h", "cpd00067"},
    "pi": {"pi", "cpd00009"}, "nadh": {"nadh", "cpd00004"},
    "nad": {"nad", "cpd00003"},
}
NAMES = {
    "atp": {"atp", "adenosine triphosphate"},
    "adp": {"adp", "adenosine diphosphate"},
    "h2o": {"h2o", "water", "h2o h2o"}, "h": {"h+", "proton"},
    "pi": {"phosphate", "orthophosphate"},
    "nadh": {"nadh", "nicotinamide adenine dinucleotide - reduced"},
    "nad": {"nad", "nad+", "nad(+)", "nicotinamide adenine dinucleotide"},
    "nadph": {"nadph", "nicotinamide adenine dinucleotide phosphate - reduced"},
    "nadp": {"nadp", "nadp+", "nadp(+)"},
}
ENERGY = [(x.upper() + "TP", {x + "tp": -1, "h2o": -1,
                            x + "dp": 1, "h": 1, "pi": 1})
          for x in ("a", "c", "g", "u", "i")]
REDOX = [(label, {red: -1, ox: 1, "h": protons}) for label, red, ox, protons in (
    ("NADH", "nadh", "nad", 1), ("NADPH", "nadph", "nadp", 1),
    ("FADH2", "fadh2", "fad", 2), ("FMNH2", "fmnh2", "fmn", 2),
    ("Q8H2", "q8h2", "q8", 2), ("MQL8", "mql8", "mqn8", 2),
    ("DMMQL8", "dmmql8", "dmmq8", 2))]


def _base_id(met):
    mid = met.id.lower()
    comp = (met.compartment or "").lower()
    # SBML may retain compartment C_c while the metabolite ID is atp_c.
    compartments = [comp, comp[2:]] if comp.startswith("c_") else [comp]
    for alias in compartments:
        for suffix in ("_" + alias, "[" + alias + "]", "@" + alias):
            if alias and mid.endswith(suffix):
                return mid[:-len(suffix)]
    return mid


def discover_targets(model):
    """Return explicit coverage; missing/ambiguous targets are never zero tests."""
    cytosol_aliases = {"c", "c0", "c_c", "c_c0", "cytosol", "cytoplasm"}
    compartments = sorted({m.compartment for m in model.metabolites
                           if (m.compartment or "").lower() in cytosol_aliases
                           or model.compartments.get(m.compartment, "").lower()
                           in {"cytosol", "cytoplasm"}})
    targets, coverage = [], []
    for comp in compartments:
        pool = [m for m in model.metabolites if m.compartment == comp]
        for label, components in ENERGY + REDOX:
            resolved, issues = {}, {}
            for key, coefficient in components.items():
                hits = [m for m in pool if _base_id(m) in ALIASES.get(key, {key})]
                if not hits:
                    hits = [m for m in pool if (m.name or "").lower() in NAMES.get(key, {key})]
                if len(hits) != 1:
                    issues[key] = "missing" if not hits else "ambiguous"
                else:
                    resolved[hits[0].id] = coefficient
            record = {"target": label + "@" + comp, "carrier": label,
                      "group": "energy" if label in {x[0] for x in ENERGY} else "redox",
                      "compartment": comp, "status": "unavailable" if issues else "available",
                      "issues": issues, "stoichiometry": resolved}
            coverage.append(record)
            if not issues:
                targets.append(record)
    if not compartments:
        coverage.append({"status": "unsupported_cytosolic_compartment"})
    return targets, coverage


def close_inputs(model):
    """Temporary zero-input conditions, independent of exchange naming/sign.

    Also relax forced nonzero fluxes in this diagnostic context so maintenance
    cannot make a zero-input test infeasible. Caller must use a model context.
    """
    for reaction in model.reactions:
        lb, ub = reaction.bounds
        if len(reaction.metabolites) == 1:
            coefficient = next(iter(reaction.metabolites.values()))
            if coefficient < 0:
                lb, ub = max(0, lb), max(0, ub)
            else:
                lb, ub = min(0, lb), min(0, ub)
        reaction.bounds = min(0, lb), max(0, ub)


def _optimize(model):
    try:
        solution = model.optimize()
        value = float(solution.objective_value) if solution.status == "optimal" else None
        if value is not None and not math.isfinite(value):
            return None, "nonfinite_objective", None
        return value, solution.status, solution
    except Exception as error:
        return None, "solver_error:" + type(error).__name__, None


def diagnose(model, target, tolerance=1e-7, witness=False):
    with model:
        close_inputs(model)
        rid = "MQC_POSTGAP_DISSIPATION"
        while rid in model.reactions:
            rid += "_"
        drain = Reaction(rid, lower_bound=0, upper_bound=1000)
        drain.add_metabolites({model.metabolites.get_by_id(mid): c
                               for mid, c in target["stoichiometry"].items()})
        model.add_reactions([drain])
        model.objective = drain
        model.objective_direction = "max"
        value, status, solution = _optimize(model)
        record = {"target": target["target"], "group": target["group"],
                  "value": value, "solver_status": status,
                  "passed": status == "optimal" and value is not None and abs(value) <= tolerance}
        if witness and value is not None and value > tolerance:
            try:
                fluxes = pfba(model).fluxes
                record["witness_method"] = "pFBA"
            except Exception:
                fluxes = solution.fluxes
                record["witness_method"] = "FBA_fallback"
            record["witness"] = [{"reaction_id": r.id, "flux": float(fluxes[r.id]),
                                  "equation": r.reaction, "bounds": list(r.bounds)}
                                 for r in model.reactions
                                 if r.id != rid and abs(fluxes[r.id]) > tolerance]
        return record


def _growth(model, biomass_id):
    with model:
        model.objective = model.reactions.get_by_id(biomass_id)
        model.objective_direction = "max"
        value, status, _ = _optimize(model)
        return {"value": value, "solver_status": status}


def _growth_viable(before, after, tolerance):
    """Reject loss of growth, not loss of a cycle-inflated growth optimum."""
    value = after.get("value")
    return (after.get("solver_status") == "optimal" and value is not None
            and math.isfinite(value)
            and (value > tolerance if before["value"] > tolerance else value >= -tolerance))


def validate_and_repair(model, biomass_id, gapfill_ids, pre_gapfill_ids,
                        tolerance=1e-7, max_edits=24):
    """Transactional guard: failed validation never leaves partial repairs.

    Greedily block a witnessed direction of a *newly added* gap reaction.
    Prefer the largest reduction in the current cycle diagnostic; ties use ID.
    Parallel cycles can retain the capped objective, so strict improvement is
    not required. Every accepted step removes a witnessed feasible direction.
    """
    targets, coverage = discover_targets(model)
    eligible = sorted((set(gapfill_ids) - set(pre_gapfill_ids)) & set(model.reactions.list_attr("id")))
    report = {"schema_version": 1, "status": "failed", "tolerance": tolerance,
              "diagnostic_drain_cap": 1000, "coverage": coverage,
              "scope": "recognized cytosolic carriers; no full thermodynamic guarantee",
              "eligible_new_gap_reactions": eligible, "attempts": [], "repairs": [],
              "before": [], "after": [], "experimental_data_used": False,
              "growth_acceptance_policy": "positive_growth_not_cycle_inflated_optimum"}
    if not biomass_id or biomass_id not in model.reactions:
        report["reason"] = "missing_biomass_objective"
        return report
    growth_before = _growth(model, biomass_id)
    report["growth_before"] = growth_before
    if growth_before["solver_status"] != "optimal" or growth_before["value"] is None:
        report["reason"] = "stored_medium_growth_not_solved"
        return report
    if not any(t["carrier"] == "ATP" for t in targets):
        report["reason"] = "ATP_target_unavailable_not_a_pass"
        return report
    # Auxiliary compartments containing no ATP pool are explicitly unassayed;
    # never combine metabolite pools from different actual compartments.
    if any(t.get("carrier") == "ATP" and t["status"] != "available"
           and t["issues"].get("atp") != "missing" for t in coverage):
        report["reason"] = "ATP_target_ambiguous_or_incomplete"
        return report
    report["before"] = [diagnose(model, t, tolerance) for t in targets]
    target_by_id = {t["target"]: t for t in targets}
    proposed = []
    with model:
        checks = report["before"]
        for iteration in range(max_edits + 1):
            report["after"] = checks
            unresolved = [d for d in checks if not d["passed"]]
            if not unresolved:
                report["status"] = "repaired" if proposed else "passed"
                report["reason"] = "all_available_targets_passed"
                break
            if any(d["solver_status"] != "optimal" or d["value"] is None for d in unresolved):
                report["reason"] = "diagnostic_solver_failure"
                break
            if iteration == max_edits:
                report["reason"] = "repair_budget_exhausted"
                break
            target = target_by_id[unresolved[0]["target"]]
            diagnostic = diagnose(model, target, tolerance, witness=True)
            attempt = {"iteration": iteration + 1, "diagnostic": diagnostic, "candidates": []}
            report["attempts"].append(attempt)
            if diagnostic["solver_status"] != "optimal" or diagnostic.get("value", 0) <= tolerance:
                report["reason"] = "nonreproducible_cycle_witness"
                break
            feasible = []
            for row in sorted(diagnostic.get("witness", []), key=lambda x: x["reaction_id"]):
                rid = row["reaction_id"]
                if rid not in eligible:
                    continue
                reaction = model.reactions.get_by_id(rid)
                lb, ub = reaction.bounds
                if not lb <= 0 <= ub:
                    continue  # Never expand a forced bound to include zero.
                bounds = (lb, 0) if row["flux"] > 0 else (0, ub)
                if bounds == reaction.bounds:
                    continue
                trial = {"reaction_id": rid, "old_bounds": [lb, ub], "new_bounds": list(bounds),
                         "blocked_direction": "forward" if row["flux"] > 0 else "reverse"}
                with model:
                    reaction.bounds = bounds
                    growth = _growth(model, biomass_id)
                    trial["growth"] = growth
                    # The pre-repair optimum may be inflated by the energy
                    # cycle just witnessed. Protect feasibility of growth,
                    # not that contaminated optimum. Rates are never fitted.
                    preserved = _growth_viable(growth_before, growth, tolerance)
                    trial["growth_preserved"] = preserved
                    if preserved:
                        trial["cycle_after"] = diagnose(model, target, tolerance)
                        if trial["cycle_after"]["solver_status"] == "optimal" and trial["cycle_after"]["value"] is not None:
                            feasible.append(trial)
                attempt["candidates"].append(trial)
            if not feasible:
                report["reason"] = "no_growth_preserving_new_gap_direction"
                break
            chosen = min(feasible, key=lambda x: (round(x["cycle_after"]["value"], 7), x["reaction_id"]))
            model.reactions.get_by_id(chosen["reaction_id"]).bounds = tuple(chosen["new_bounds"])
            proposed.append(chosen)
            attempt["chosen_reaction"] = chosen["reaction_id"]
            checks = [diagnose(model, t, tolerance) for t in targets]
        report["growth_after_trial"] = _growth(model, biomass_id)
    report["proposed_repairs"] = proposed
    if report["status"] != "failed":
        for change in proposed:
            model.reactions.get_by_id(change["reaction_id"]).bounds = tuple(change["new_bounds"])
        report["repairs"] = proposed
    else:
        report["trial_changes_rolled_back"] = True
    return report


def run_postgap_control(model_info, model, model_control_info, cfg, pre_gapfill_ids):
    """Persist auditable status and prevent legacy RuntimeError success fallback."""
    try:
        environment = validate_growth_environment(model, getattr(cfg, "growth_environment", "preserve"))
        report = validate_and_repair(model, model_info.get("initial_rxn", {}).get("biomass_rxn_id"),
                                     model_info.get("gaped_rxn", []), pre_gapfill_ids)
        report["growth_environment"] = environment
    except Exception as error:
        # Trials are protected by model contexts; never swallow an unknown error
        # as an apparently clean final model.
        report = {"status": "failed", "reason": "unexpected_validation_error",
                  "error": type(error).__name__ + ": " + str(error), "repairs": []}
    model_info["postgap_validation"] = report
    model_control_info["postgap_validation"] = report
    for change in report["repairs"]:
        for entry in model_info.get("reactions", []):
            if entry.get("id") == change["reaction_id"]:
                entry["bounds"] = tuple(change["new_bounds"])
                entry["postgap_bounds_revision"] = change
    Path(cfg.output_dir).mkdir(parents=True, exist_ok=True)
    (Path(cfg.output_dir) / "postgap_validation.json").write_text(
        json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False), encoding="utf-8")
    if report["status"] == "failed":
        model_control_info["final_model_status"] = "postgap_validation_failed"
        for key in ("check_energy_production", "check_reducing_equivalents_production"):
            model_control_info.setdefault(key, {})["postgap_validation_status"] = "failed_or_incomplete"
        raise PostGapValidationError("Final post-gapfill validation failed: " + report["reason"])
    model_control_info["final_model_status"] = "postgap_validated"
    return report
