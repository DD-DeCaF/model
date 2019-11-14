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


def minimize_distance(model, biomass_reaction, growth_rate, fluxomics):
    """Replaces fluxomics measurements with the minimized distance"""
    index = []
    observations = []
    uncertainties = []

    if not growth_rate:
        raise ValueError(
            "Expected measurements to contain an objective "
            "constraint as measured growth rate"
        )

    # Trust the growth rate over the measurements. Meaning, constrain the
    # biomass reaction to the observed values instead of simply including it in
    # the measurements to be minimized.
    if growth_rate["uncertainty"]:
        lower_bound = growth_rate["measurement"] - growth_rate["uncertainty"]
        upper_bound = growth_rate["measurement"] + growth_rate["uncertainty"]
    else:
        lower_bound = growth_rate["measurement"]
        upper_bound = growth_rate["measurement"]
    model.reactions.get_by_id(biomass_reaction).bounds = (lower_bound, upper_bound)

    for measure in fluxomics:
        index.append(measure["identifier"])
        observations.append(measure["measurement"])
        # TODO: How to implement uncertainty here?
        uncertainties.append(1)

    observations = pd.Series(index=index, data=observations)
    uncertainties = pd.Series(index=index, data=uncertainties)

    solution = adjust_fluxes2model(model, observations, uncertainties)
    for reaction, minimized_distance in solution.fluxes.iteritems():
        if reaction == biomass_reaction:
            growth_rate["measurement"] = minimized_distance
        for measure in fluxomics:
            if reaction == measure.get("identifier"):
                measure["measurement"] = minimized_distance
                measure["uncertainty"] = 0  # TODO: Confirm that this is correct
    return growth_rate, fluxomics


def adjust_fluxes2model(
    model, observations, uncertainties=None, linear=True, big_m=1e05
):
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
    data.loc[
        data[weight_col].isnull()
        | np.isinf(data[weight_col])
        | (data[weight_col] == 0),
        weight_col,
    ] = 1
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
                    ub=0,
                    name="forward_pos_" + rxn_id,
                )
                forward_neg = prob.Constraint(
                    rxn.flux_expression - flux - big_m * (1 - direction) - dist,
                    ub=0,
                    name="forward_neg_" + rxn_id,
                )
                reverse_pos = prob.Constraint(
                    (-flux) - rxn.flux_expression - big_m * direction - dist,
                    ub=0,
                    name="reverse_pos_" + rxn_id,
                )
                reverse_neg = prob.Constraint(
                    rxn.flux_expression - (-flux) - big_m * direction - dist,
                    ub=0,
                    name="reverse_neg_" + rxn_id,
                )
                if linear:
                    new_obj += dist / weight
                else:
                    new_obj += (dist / weight) ** 2
                to_add.extend(
                    [
                        direction,
                        dist,
                        forward_pos,
                        forward_neg,
                        reverse_pos,
                        reverse_neg,
                    ]
                )
            except KeyError:
                logger.warning(
                    f"Reaction '{rxn_id}' not found in the model. " f"Ignored."
                )
        model.add_cons_vars(to_add)
        model.objective = prob.Objective(new_obj, direction="min")
        solution = model.optimize(raise_error=True)
    return solution


def flexibilize_proteomics(model, biomass_reaction, growth_rate, proteomics):
    """
    Replace proteomics measurements with a set that enables the model to grow. Proteins
    are removed from the set iteratively based on sensitivity analysis (shadow prices).

    Parameters
    ----------
    model: cobra.Model
        The enzyme-constrained model.
    biomass_reaction: str
        The id of the biomass reaction in the given model.
    proteomics: list(dict)
        List of measurements matching the `Proteomics` schema.
    growth_rate: dict
        Growth rate, matching the `GrowthRate` schema.

    Returns
    -------
    growth_rate: dict
        New growth rate (will change if the model couldn't grow at the inputted value).
    proteomics: list(dict)
        Filtered list of proteomics.
    """

    # reset growth rate in model:
    model.reactions.get_by_id(biomass_reaction).bounds = (0, 1000)

    # compute measurements to constrain with:
    measurements = pd.Series()
    for protein in proteomics:
        protein_id = protein["identifier"]
        value = protein["measurement"] + protein["uncertainty"]
        measurements = measurements.append(pd.Series(data=[value], index=[protein_id]))

    # constrain the model with all proteins and optimize:
    limit_proteins(model, measurements)
    solution = model.optimize()
    new_growth_rate = solution.objective_value

    # while the model cannot grow to the desired level, remove the protein with
    # the higher shadow price:
    desired_growth = growth_rate["measurement"]
    prots_to_remove = []
    while new_growth_rate < desired_growth and not measurements.empty:
        # get most influential protein in model:
        top_protein = top_protein_shadow_prices(solution, measurements.index)
        top_protein = top_protein.index[0]

        # update data: append protein to list, remove from current dataset and
        # open the corresponding upper bound:
        prots_to_remove.append(top_protein)
        measurements.pop(top_protein)
        rxn = model.reactions.get_by_id("prot_{}_exchange".format(top_protein))
        rxn.bounds = (0, 1000)

        # rinse and repeat:
        limit_proteins(model, measurements)
        solution = model.optimize()
        new_growth_rate = solution.objective_value

    # update growth rate if optimization was not succesful:
    if new_growth_rate < desired_growth:
        growth_rate["measurement"] = new_growth_rate

    # update proteomics by removing flexibilized proteins:
    for protein in prots_to_remove:
        index = next((index for (index, d) in enumerate(proteomics) if d["identifier"] == protein), None)
        proteomics.pop(index)

    return growth_rate, proteomics


def limit_proteins(model, measurements):
    """Apply proteomics measurements to model.

    Parameters
    ----------
    model: cobra.Model
        The enzyme-constrained model.
    measurements : pd.Series
        Protein abundances in mmol / gDW.
    """
    for protein_id, measure in measurements.iteritems():
        try:
            rxn = model.reactions.get_by_id("prot_{}_exchange".format(protein_id))
        except KeyError:
            pass
        else:
            # update only upper_bound (as enzymes can be unsaturated):
            rxn.bounds = (0, measure)
    return


def top_protein_shadow_prices(solution, proteins, top=1):
    """
    Retrieves shadow prices of proteins from solution and ranks them from most
    to least sensitive in the model.

    Parameters
    ----------
    solution: cobra.Solution
        The usual Solution object returned by model.optimize().
    proteins: iterable of strings
        Protein IDs of the proteins in measurements.
    top: int
        The number of proteins to be returned.

    Returns
    -------
    shadow_proteins: pd.Series
        Ranked shadow prices.
    """
    # TODO: refactor for speed
    shadow_proteins = pd.Series()
    for protein in proteins:
        for key, value in solution.shadow_prices.items():
            if protein in key:
                shadow_proteins[protein] = value
    return shadow_proteins.sort_values()[:top]
