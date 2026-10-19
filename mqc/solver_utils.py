"""Bounded numerical recovery; never discard constraints or accept bad fluxes."""
import logging
import math

from cobra.exceptions import OptimizationError
from cobra.flux_analysis import pfba as cobra_pfba
from cobra.flux_analysis.parsimonious import add_pfba
from cobra.core.solution import get_solution
from cobra.util.solver import linear_reaction_coefficients

LOG = logging.getLogger(__name__)


def exceeds_yield_limit(actual, limit, atol=1e-7, rtol=1e-6):
    """Compare unrounded fluxes; nonfinite results are errors, not a pass."""
    actual, limit = float(actual), float(limit)
    if not math.isfinite(actual) or not math.isfinite(limit):
        raise OptimizationError('Nonfinite metabolite yield or theoretical limit')
    return actual - limit > atol + rtol * max(abs(actual), abs(limit))


def cold_copy(model):
    """Rebuild the complete optlang problem, including custom constraints."""
    model.solver.update()
    data = model.solver.to_json()
    result = model.copy()
    result.solver = model.problem.Model.from_json(data)
    old, new = model.solver.configuration, result.solver.configuration
    for key in ('timeout', 'presolve', 'verbosity'):
        setattr(new, key, getattr(old, key))
    for key in ('feasibility', 'optimality', 'integrality'):
        if hasattr(old.tolerances, key) and hasattr(new.tolerances, key):
            setattr(new.tolerances, key, getattr(old.tolerances, key))
    # This recovery is used in per-model workers, not a new parallel layer.
    if hasattr(result.solver.problem, 'parameters'):
        result.solver.problem.parameters.threads.set(1)
    return result


def validate_primals(model, tolerance=1e-6):
    """Check all LP bounds, including extra (non-stoichiometric) constraints."""
    if model.solver.status != 'optimal':
        raise OptimizationError('Recovery solve is not optimal: ' + model.solver.status)
    worst = 0.0
    for items, values in [(model.solver.variables, model.solver.primal_values),
                          (model.solver.constraints, model.solver.constraint_values)]:
        for item in items:
            value = values[item.name]
            if value is None or not math.isfinite(value):
                raise OptimizationError('Nonfinite primal: ' + item.name)
            worst = max(worst, 0 if item.lb is None else item.lb - value,
                        0 if item.ub is None else value - item.ub)
    if worst > tolerance:
        raise OptimizationError('Recovery violates an LP bound by %.9g' % worst)
    return worst


def _checked_pfba(model, fraction):
    # CPLEX invalidates its primal state when the pFBA context is removed.
    # Validate before removing the temporary objective/constraint, not after.
    with model:
        add_pfba(model, fraction_of_optimum=fraction)
        model.slim_optimize(error_value=None)
        residual = validate_primals(model)
        return get_solution(model), residual


def solve_with_cold_retry(model):
    """Return working model, value and provenance; retry known LP failures once.

    A second failure is propagated. Optional callers may explicitly skip their
    extension, but must not label an unresolved optimization a biological pass.
    """
    try:
        value = model.slim_optimize(error_value=None)
        if not math.isfinite(value):
            raise OptimizationError('Nonfinite primary objective')
        return model, value, {'method': 'native', 'status': 'optimal'}
    except OptimizationError as error:
        recovered = cold_copy(model)
        value = recovered.slim_optimize(error_value=None)
        residual = validate_primals(recovered)
        if not math.isfinite(value):
            raise OptimizationError('Nonfinite recovered objective')
        info = dict(method='complete_lp_cold_retry', status='optimal',
                    initial_error=repr(error), max_bound_violation=residual,
                    variables=len(recovered.variables), constraints=len(recovered.constraints))
        LOG.warning('MQC numerical recovery: %s', info)
        return recovered, value, info


def robust_pfba(model, fraction_of_optimum=1.0, objective=None, reactions=None):
    """Preserve ordinary pFBA; recover only an OptimizationError on a copy.

    Try exact cold-start first. For a positive maximizing objective only, a
    final 1e-6 relative objective relaxation is permitted. Check that objective
    error, Sv, reaction bounds AND custom LP constraints satisfy tolerances.
    No model bounds, objective or solver on the caller's model are replaced.
    """
    try:
        return cobra_pfba(model, fraction_of_optimum, objective, reactions)
    except OptimizationError as original:
        if fraction_of_optimum != 1.0 or objective is not None or reactions is not None:
            raise
        work, optimum, primary_info = solve_with_cold_retry(cold_copy(model))
        coefficients = {r.id: c for r, c in linear_reaction_coefficients(work).items()}
        if not coefficients:
            raise original
        fraction = 1.0
        try:
            solution, residual = _checked_pfba(work, fraction)
        except OptimizationError:
            if work.objective_direction != 'max' or optimum <= 0:
                raise
            fraction = 1.0 - 1e-6
            solution, residual = _checked_pfba(work, fraction)
        primary = sum(c * float(solution.fluxes[r]) for r, c in coefficients.items())
        if not math.isfinite(primary) or abs(primary - optimum) > 1e-9 + 1e-6 * abs(optimum):
            raise OptimizationError('pFBA recovery did not preserve the primary objective')
        sv = {m.id: 0.0 for m in work.metabolites}
        bound_error = 0.0
        for reaction in work.reactions:
            value = float(solution.fluxes[reaction.id])
            if not math.isfinite(value):
                raise OptimizationError('Nonfinite recovered pFBA flux')
            bound_error = max(bound_error, reaction.lower_bound - value, value - reaction.upper_bound)
            for metabolite, coefficient in reaction.metabolites.items():
                sv[metabolite.id] += coefficient * value
        if max(bound_error, max(map(abs, sv.values()), default=0)) > 1e-6:
            raise OptimizationError('pFBA recovery failed flux feasibility checks')
        LOG.warning('MQC pFBA recovery: original=%r fraction=%s optimum=%s actual=%s '
                    'max_LP_violation=%s primary=%s', original, fraction, optimum,
                    primary, residual, primary_info)
        return solution
