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
import re

import numpy as np
import pandas as pd
from optlang.symbolics import Zero
from math import isnan


logger = logging.getLogger(__name__)
prot_pat = re.compile(r"prot_([A-Za-z0-9]+)__91__c__93__")


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
    Replace proteomics measurements with the set of proteomics measures that enables
    the model to achieve growth.

    Proteins are removed from the set iteratively based on sensitivity analysis.

    Parameters
    ----------
    model: cobra.Model
    biomass_reaction: str
        The id of the biomass reaction in the goperationsiven model.
    proteomics: list(dict)
        List of measurements matching the `Proteomics` schema.
    growth_rate: dict
        Growth rate, matching the `GrowthRate` schema.

    Returns
    -------
    growth_rate: dict
    proteomics: list(dict)
    """
    # TODO: this whole thing about growth rate could be refactored, since it's exactly
    # the same as in `minimize_distance`
    if not growth_rate:
        raise ValueError(
            "Expected measurements to contain an objective "
            "constraint as measured growth rate"
        )

    if growth_rate["uncertainty"]:
        lower_bound = growth_rate["measurement"] - growth_rate["uncertainty"]
        upper_bound = growth_rate["measurement"] + growth_rate["uncertainty"]
    else:
        lower_bound = growth_rate["measurement"]
        upper_bound = growth_rate["measurement"]

    model.reactions.get_by_id(biomass_reaction).bounds = (lower_bound, upper_bound)

    index, observations = [], []

    for measure in proteomics:
        index.append(measure["identifier"])
        # in protein data, upper_bound is the unique constraint, so it contains the uncertainty
        observations.append(measure["measurement"] + measure["uncertainty"])

    observations = pd.Series(index=index, data=observations)

    solution, new_growth_rate = ensure_proteomics_growing(model, observations)
    if new_growth_rate:
        growth_rate["measurement"] = new_growth_rate
    for reaction, new_measurement in solution.iteritems():
        for measure in proteomics:
            if reaction == measure.get("identifier"):
                measure["measurement"] = new_measurement
                measure["uncertainty"] = 0  # TODO: Confirm that this is correct
    return growth_rate, proteomics


def limit_proteins(model, measured_mmolgdw):
    """Apply proteomics measurements to model.

    Apply measurements in the form of measured proteins (mmol/gDW).

    Parameters
    ----------
    model: cobra.Model
    measured_mmolgdw : pd.Series
        Protein abundances in mmol / gDW

    References
    ----------
    .. [1] Benjamin J. Sanchez, Cheng Zhang, Avlant Nilsson, Petri-Jaan Lahtvee, Eduard J. Kerkhoven, Jens Nielsen (
       2017). Improving the phenotype predictions of a yeast genome-scale metabolic model by incorporating enzymatic
       constraints. [Molecular Systems Biology, 13(8): 935, http://www.dx.doi.org/10.15252/msb.20167411
    """
    # update upper_bound
    for protein_id, measure in measured_mmolgdw.iteritems():
        try:
            rxn = model.reactions.get_by_id("prot_{}_exchange".format(protein_id))
            rxn.annotation["uniprot"] = protein_id
        except KeyError:
            pass
        else:
            rxn.bounds = 0, measure
    # flexibilize proteomics so model grows
    return


def ensure_proteomics_growing(model, measured_mmolgdw):
    """
    Optimize the model and flexiblize proteins until it grows

    Parameters
    ----------
    model: cobra.Model
    measured_mmolgdw: list(dict)
        List of measurements matching the `Proteomics` schema.

    Returns
    -------
    measured_mmolgdw: list(dict)
    new_growth_rate: float. False if it wasn't changed
    """
    # first, get shadow prices of the unconstrained model
    solution = model.optimize()
    ordered_proteins = list(
        top_protein_shadow_prices(solution, measured_mmolgdw.index, -1).index
    )
    # second, constrain the model
    with model as mod:
        limit_proteins(mod, measured_mmolgdw)
        obj_value = mod.slim_optimize()
    tolerance = 1e-7
    obj_value = 0 if math.isnan(obj_value) else obj_value
    new_growth_rate = False

    # while the model can't grow, unconstrain the protein with the higher shadow price
    while obj_value < tolerance and not measured_mmolgdw.empty:
        uniprot_id = re.sub(prot_pat, r"\1", ordered_proteins[0])
        if uniprot_id in measured_mmolgdw:
            del measured_mmolgdw[uniprot_id]
        ordered_proteins.pop(0)
        with model as mod:
            # this can be change to simply affect the upper_bound
            # but this way is more clear
            limit_proteins(mod, measured_mmolgdw)
            obj_value = mod.slim_optimize()
        new_growth_rate = obj_value = 0 if math.isnan(obj_value) else obj_value

    return measured_mmolgdw, new_growth_rate


def top_protein_shadow_prices(model_solution, set_proteins, top=1):
    """
    Retrieves `top` of proteins in `set_proteins` in terms of influence in the objective
    function (shadow prices).

    Parameters
    ----------
    model_solution: cobra.Solution
        the usual Solution object returned by model.optimize()
    set_proteins: iterable of strings
        Uniprot IDs of the proteins in measurements
    top: int
        the number of proteins to be returned
    """
    set_as_metabolites = {"prot_{}__91__c__93__".format(prot) for prot in set_proteins}
    shadow_pr = model_solution.shadow_prices
    return shadow_pr.loc[shadow_pr.index.isin(set_as_metabolites)].sort_values()[:top]
