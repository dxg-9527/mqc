"""Explicit growth-environment constraints, independent of experimental rates.

Anaerobic QC must gap-fill with the oxygen/ROS uptake restrictions in place,
not discover only after export that the selected pathway needs oxygen. The
default preserves legacy behavior; this policy is never inferred from an ID.
"""


class GrowthEnvironmentError(ValueError):
    """A declared physiological condition is unsupported or was violated."""


def input_closed_bounds(reaction):
    """Block input of a single-metabolite boundary; retain allowed secretion.

    The convention comes from the stoichiometric sign, not the reaction ID.
    A forced input becomes zero; previously absent secretion is not invented.
    """
    if len(reaction.metabolites) != 1:
        raise ValueError("Input closure requires a single-metabolite boundary")
    coefficient = next(iter(reaction.metabolites.values()))
    lb, ub = reaction.bounds
    return (max(0, lb), max(0, ub)) if coefficient < 0 else (min(0, lb), min(0, ub))


def _oxygen_or_ros(met):
    # O2 includes molecular oxygen and the superoxide formula O2 (charge -1).
    # Water, sulfate and other oxygen-containing nutrients are not excluded.
    if met.formula in {"O2", "H2O2"}:
        return True
    names = {"oxygen", "o2", "o2 o2", "hydrogen peroxide", "h2o2", "superoxide",
             "superoxide anion", "o2-", "o2s"}
    if (met.name or "").strip().lower() in names:
        return True
    mid = met.id.lower()
    return any(mid == base or mid.startswith(base + "_") or mid.startswith(base + "[")
               for base in ("o2", "o2s", "h2o2", "cpd00007", "cpd00025", "cpd00532"))


def oxygen_boundaries(model):
    return [r for r in model.reactions if len(r.metabolites) == 1
            and _oxygen_or_ros(next(iter(r.metabolites)))]


def apply_growth_environment(model, mode="preserve"):
    if mode not in {"preserve", "anaerobic"}:
        raise GrowthEnvironmentError("growth_environment must be preserve or anaerobic")
    report = {"mode": mode, "changes": [], "experimental_rates_used": False}
    if mode == "preserve":
        return report
    for reaction in oxygen_boundaries(model):
        before = reaction.bounds
        after = input_closed_bounds(reaction)
        if after != before:
            reaction.bounds = after
            report["changes"].append({"reaction_id": reaction.id, "before": list(before),
                                      "after": list(after), "equation": reaction.reaction})
    return report


def validate_growth_environment(model, mode="preserve"):
    if mode not in {"preserve", "anaerobic"}:
        raise GrowthEnvironmentError("growth_environment must be preserve or anaerobic")
    violations = [r.id for r in oxygen_boundaries(model) if r.bounds != input_closed_bounds(r)] if mode == "anaerobic" else []
    if violations:
        raise GrowthEnvironmentError("Declared anaerobic uptake restrictions were reopened: " + ", ".join(violations))
    return {"mode": mode, "status": "passed", "checked_boundary_ids":
            [r.id for r in oxygen_boundaries(model)] if mode == "anaerobic" else []}
