import logging

from numpy import isinf
from optlang.symbolics import Zero
from pandas import Series


LOGGER = logging.getLogger(__name__)


"""This module may one day be replaced with http://driven.bio/"""


def adjust_fluxes2model(model, observations, uncertainties=None, linear=True,
                        big_m=1E05):
    """
    Minimize the distance to observed fluxes accounting for multiple directions.

    If your observations include uncertainties the objective function, i.e.,
    minimizing the distance to the observations, is weighted by the inverse
    of the uncertainties.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model
    observations : pandas.Series
        The observed fluxes. The index should contain reaction identifiers.
    uncertainties : pandas.Series, optional
        The uncertainties of the individual measurements, e.g., standard
        error. The index of the series should correspond at least partially
        to the ``observations``.
    linear : bool, optional
        Whether to minimize the linear or quadratic distance.
    big_m : float, optional
        Big M method value. This is used to resolve greater than inequalities
        and should be an adequately large number.

    Returns
    -------
    cobra.Solution

    """
    flux_col = "flux"
    weight_col = "weight"
    if uncertainties is None:
        data = observations.to_frame().join(Series([], name=weight_col))
    else:
        uncertainties.name = weight_col
        data = observations.to_frame().join(uncertainties)
    data.columns = [flux_col, weight_col]
    # replace missing and zero values
    data.loc[data[weight_col].isnull() | isinf(data[weight_col]) |
             (data[weight_col] == 0), weight_col] = 1
    prob = model.problem
    to_add = list()
    new_obj = Zero
    with model:
        for rxn_id, flux, weight in data[[flux_col, weight_col]].itertuples():
            try:
                rxn = model.reactions.get_by_id(rxn_id)
                direction = prob.Variable("direction_" + rxn_id, type="binary")
                dist = prob.Variable("dist_" + rxn_id)
                forward_pos = prob.Constraint(
                    flux - rxn.flux_expression - big_m * (1 - direction) - dist,
                    ub=0, name="forward_pos_" + rxn_id)
                forward_neg = prob.Constraint(
                    rxn.flux_expression - flux - big_m * (1 - direction) - dist,
                    ub=0, name="forward_neg_" + rxn_id)
                reverse_pos = prob.Constraint(
                    (-flux) - rxn.flux_expression - big_m * direction - dist,
                    ub=0, name="reverse_pos_" + rxn_id)
                reverse_neg = prob.Constraint(
                    rxn.flux_expression - (-flux) - big_m * direction - dist,
                    ub=0, name="reverse_neg_" + rxn_id)
                if linear:
                    new_obj += (dist / weight)
                else:
                    new_obj += (dist / weight) ** 2
                to_add.extend([direction, dist, forward_pos, forward_neg,
                               reverse_pos, reverse_neg])
            except KeyError:
                LOGGER.warning(f"Reaction '{rxn_id}' not found in the model. "
                               f"Ignored.")
        model.add_cons_vars(to_add)
        model.objective = prob.Objective(new_obj, direction="min")
        solution = model.optimize(raise_error=True)
    return solution