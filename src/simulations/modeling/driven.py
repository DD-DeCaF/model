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
    lb, ub = bounds(growth_rate["measurement"], growth_rate["uncertainty"])
    model.reactions.get_by_id(biomass_reaction).bounds = (lb, ub)

    for measure in fluxomics:
        index.append(measure["identifier"])
        observations.append(measure["measurement"])
        # TODO: How to implement uncertainty here?
        uncertainties.append(1)

    observations = pd.Series(index=index, data=observations)
    uncertainties = pd.Series(index=index, data=uncertainties)

    solution = adjust_fluxes2model(model, observations, uncertainties)
    for reaction, minimized_distance in solution.fluxes.items():
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
    growth_rate: dict
        Growth rate, matching the `GrowthRate` schema.
    proteomics: list(dict)
        List of measurements matching the `Proteomics` schema.

    Returns
    -------
    growth_rate: dict
        New growth rate (will change if the model couldn't grow at the inputted value).
    proteomics: list(dict)
        Filtered list of proteomics.
    """

    # reset growth rate in model:
    model.reactions.get_by_id(biomass_reaction).bounds = (0, 1000)

    # build a table with protein ids, met ids in model and values to constrain with:
    prot_df = pd.DataFrame()
    for protein in proteomics:
        protein_id = protein["identifier"]
        lb, ub = bounds(protein["measurement"], protein["uncertainty"])
        for met in model.metabolites.query(lambda m: protein_id in m.id):
            new_row = pd.DataFrame(
                data={"met_id": met.id, "value": ub}, index=[protein_id]
            )
            prot_df = prot_df.append(new_row)

    # constrain the model with all proteins and optimize:
    limit_proteins(model, prot_df["value"])
    solution = model.optimize()
    new_growth_rate = solution.objective_value

    # while the model cannot grow to the desired level, remove the protein with
    # the highest shadow price:
    minimal_growth, ub = bounds(growth_rate["measurement"], growth_rate["uncertainty"])
    prots_to_remove = []
    while new_growth_rate < minimal_growth and not prot_df.empty:
        # get most influential protein in model:
        top_protein = top_shadow_prices(solution, list(prot_df["met_id"]))
        top_protein = top_protein.index[0]
        top_protein = prot_df.index[prot_df["met_id"] == top_protein][0]

        # update data: append protein to list, remove from current dataframe and
        # increase the corresponding upper bound to +1000:
        prots_to_remove.append(top_protein)
        prot_df = prot_df.drop(labels=top_protein)
        limit_proteins(model, pd.Series(data=[1000], index=[top_protein]))

        # re-compute solution:
        solution = model.optimize()
        if solution.objective_value == new_growth_rate:  # the algorithm is stuck
            break
        new_growth_rate = solution.objective_value

    # update growth rate if optimization was not successful:
    if new_growth_rate < minimal_growth:
        if growth_rate["uncertainty"]:
            growth_rate["measurement"] = new_growth_rate + growth_rate["uncertainty"]
        else:
            growth_rate["measurement"] = new_growth_rate

    # update proteomics by removing flexibilized proteins:
    for protein in prots_to_remove:
        index = next(
            (
                index
                for (index, dic) in enumerate(proteomics)
                if dic["identifier"] == protein
            ),
            None,
        )
        del proteomics[index]

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
    for protein_id, measure in measurements.items():
        try:
            rxn = model.reactions.get_by_id(f"prot_{protein_id}_exchange")
        except KeyError:
            pass
        else:
            # update only upper_bound (as enzymes can be unsaturated):
            rxn.bounds = (0, measure)
    return


def top_shadow_prices(solution, met_ids, top=1):
    """
    Retrieves shadow prices for a list of metabolites from the solution and ranks
    them from most to least sensitive in the model.

    Parameters
    ----------
    solution: cobra.Solution
        The usual Solution object returned by model.optimize().
    met_ids: iterable of strings
        Subset of metabolite IDs from the model.
    top: int
        The number of metabolites to be returned.

    Returns
    -------
    shadow_pr: pd.Series
        Top shadow prices, ranked.
    """
    shadow_pr = solution.shadow_prices
    shadow_pr = shadow_pr.loc[shadow_pr.index.isin(met_ids)]
    return shadow_pr.sort_values()[:top]


def bounds(measurement, uncertainty):
    """Return resolved bounds based on measurement and uncertainty"""
    if uncertainty:
        return measurement - uncertainty, measurement + uncertainty
    else:
        return measurement, measurement
