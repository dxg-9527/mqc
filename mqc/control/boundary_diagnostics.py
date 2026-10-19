"""Direction-aware boundary handling for temporary diagnostic conditions."""
from contextlib import contextmanager
import math
import re

from mqc.control.growth_environment import input_closed_bounds


def unresolved_formula(metabolite):
    formula = metabolite.formula
    if formula is None:
        return True
    value = str(formula).strip()
    return (value.lower() in {"", "nan", "none", "null"}
            or any(element in {"R", "X"} for element in re.findall(r"[A-Z][a-z]?", value)))


def boundary_input_flux(reaction, flux):
    """Positive means material enters the model, irrespective of ID or sign."""
    if len(reaction.metabolites) != 1:
        raise ValueError("Expected a single-metabolite boundary")
    return next(iter(reaction.metabolites.values())) * float(flux)


def unresolved_boundary_source(reaction, flux, tolerance=1e-6):
    return (len(reaction.metabolites) == 1
            and math.isfinite(float(flux))
            and boundary_input_flux(reaction, flux) > tolerance
            and unresolved_formula(next(iter(reaction.metabolites))))


@contextmanager
def without_unresolved_inputs(model):
    """Temporarily remove unclassified nutrient inputs, not their secretion.

    Used only for the zero-carbon net-production assay. Restore these boundary
    bounds even on failure; retain unrelated internal-reaction corrections.
    """
    changes = []
    try:
        for reaction in model.reactions:
            if len(reaction.metabolites) != 1 or not unresolved_formula(next(iter(reaction.metabolites))):
                continue
            old, new = reaction.bounds, input_closed_bounds(reaction)
            if old != new:
                changes.append({"reaction_id": reaction.id, "old_bounds": list(old),
                                "diagnostic_bounds": list(new)})
                reaction.bounds = new
        yield changes
    finally:
        for change in changes:
            model.reactions.get_by_id(change["reaction_id"]).bounds = tuple(change["old_bounds"])
