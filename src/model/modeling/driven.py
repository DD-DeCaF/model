# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging

import numpy as np
import pandas as pd
from optlang.symbolics import Zero


logger = logging.getLogger(__name__)


"""This module may one day be replaced with http://driven.bio/"""


def minimize_distance(model, biomass_reaction, growth_rate, measurements):
    """Replaces fluxomics measurements with the minimized distance"""
    index = []
    observations = []
    uncertainties = []

    if not growth_rate:
        raise ValueError("Expected measurements to contain an objective "
                         "constraint as measured growth rate")

    # Trust the growth rate over the measurements. Meaning, constrain the
    # biomass reaction to the observed values instead of simply including it in
    # the measurements to be minimized.
    # TODO (Ali Kaafarani): Support for uncertainties or multiple observations
    if len(growth_rate['measurements']) != 1:
        raise NotImplementedError("Cannot handle multiple growth rate measurements yet")
    model.reactions.get_by_id(biomass_reaction).bounds = (growth_rate['measurements'][0], growth_rate['measurements'][0])

    for measure in [m for m in measurements if m['type'] == 'reaction']:
        index.append(measure['id'])
        observations.append(np.nanmean(measure['measurements']))
        if len(measure['measurements']) >= 3:
            uncertainties.append(np.nanstd(measure['measurements'], ddof=1))
        else:
            uncertainties.append(1)

    observations = pd.Series(index=index, data=observations)
    uncertainties = pd.Series(index=index, data=uncertainties)

    solution = adjust_fluxes2model(model, observations, uncertainties)
    for reaction, minimized_distance in solution.fluxes.iteritems():
        if reaction == biomass_reaction:
            growth_rate['measurements'] = [minimized_distance]
        for measurement in measurements:
            if reaction == measurement.get('id'):
                measurement['measurements'] = [minimized_distance]
    return growth_rate, measurements


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
        data = observations.to_frame().join(pd.Series([], name=weight_col))
    else:
        uncertainties.name = weight_col
        data = observations.to_frame().join(uncertainties)
    data.columns = [flux_col, weight_col]
    # replace missing and zero values
    data.loc[data[weight_col].isnull() | np.isinf(data[weight_col]) |
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
                logger.warning(f"Reaction '{rxn_id}' not found in the model. "
                               f"Ignored.")
        model.add_cons_vars(to_add)
        model.objective = prob.Objective(new_obj, direction="min")
        solution = model.optimize(raise_error=True)
    return solution
