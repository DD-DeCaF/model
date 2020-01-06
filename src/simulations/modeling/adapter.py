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

import json
import logging
from collections import namedtuple

from cobra import Configuration, Metabolite, Reaction
from cobra.io.dict import reaction_to_dict
from cobra.medium.boundary_types import find_external_compartment

from simulations.exceptions import CompartmentNotFound, MetaboliteNotFound, PartNotFound
from simulations.ice_client import ICE
from simulations.modeling.cobra_helpers import (
    find_metabolite,
    get_exchange_reaction,
    parse_bigg_compartment,
)
from simulations.modeling.driven import flexibilize_proteomics, minimize_distance
from simulations.modeling.gnomic_helpers import feature_id


logger = logging.getLogger(__name__)
ice = ICE()

with open("data/salts.json") as file_:
    SALTS = json.load(file_)


def apply_medium(model, is_ec_model, medium):
    """
    Apply a medium to a metabolic model.

    The medium is applied by adding salt mappings, resolving the corresponding
    exchange reaction for the resulting medium compounds, setting their uptake
    rate, and then applying that to the model.

    Parameters
    ----------
    model: cobra.Model
    is_ec_model: bool
        A boolean indicating if the model is enzyme-constrained.
    medium: list(dict)
        The medium definition, a list of dicts matching the `MediumCompound`
        schema.

    Returns
    -------
    tuple (operations, warnings, errors)
        Operations is a list of model operations necessary to apply the medium
        to the model. See also the `Operations` schema.
        Warnings is a list of human-readable strings of potential issues.
        If errors is not an empty list, it was not possible to apply the medium.
        Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    warnings = []
    errors = []

    # Convert the list of dicts to a set of namedtuples to avoid duplicates, as
    # looking up metabolites in the model is a somewhat expensive operation.
    Compound = namedtuple("Compound", ["id", "namespace"])
    medium = set(Compound(id=c["identifier"], namespace=c["namespace"]) for c in medium)

    # Detect salt compounds and split them into their ions and metals
    for compound in medium.copy():  # Make a copy to be able to mutate the original list
        if compound.id in SALTS:
            salt = SALTS[compound.id]
            logger.info(
                f"Replacing {compound.id} with ions: {salt['ions']} and metals: "
                f"{salt['metals']}"
            )
            medium.remove(compound)
            medium.update([Compound(id=ion, namespace="chebi") for ion in salt["ions"]])
            medium.update(
                [Compound(id=metal, namespace="chebi") for metal in salt["metals"]]
            )

            if salt["ions_missing_smiles"]:
                warning = (
                    f"Unable to add ions, smiles id could not be mapped: "
                    f"{salt['ions_missing_smiles']}"
                )
                warnings.append(warning)
                logger.warning(warning)
            if salt["metals_missing_inchi"]:
                warning = (
                    f"Unable to add metals; inchi string could not be mapped: "
                    f"{salt['metals_missing_inchi']}"
                )
                warnings.append(warning)
                logger.warning(warning)

    # Add trace metals
    medium.update(
        [
            Compound(id="CHEBI:25517", namespace="chebi"),
            Compound(id="CHEBI:25368", namespace="chebi"),
        ]
    )

    try:
        extracellular = find_external_compartment(model)
    except RuntimeError as error:
        # cobrapy throws RuntimeError if it for any reason is unable to find an
        # external compartment. See:
        # https://github.com/opencobra/cobrapy/blob/95d920d135fa824e6087f1fcbc88d50882da4dab/cobra/medium/boundary_types.py#L26
        message = (
            f"Cannot find an external compartment in model {model.id}: {str(error)}"
        )
        errors.append(message)
        logger.error(message)
        # Cannot continue without knowing the external compartment, so
        # immediately return the error.
        return operations, warnings, errors

    # Create a map of exchange reactions and corresponding fluxes to apply to
    # the medium.
    medium_mapping = {}
    for compound in medium:
        try:
            extracellular_metabolite = find_metabolite(
                model, compound.id, compound.namespace, extracellular
            )
        except MetaboliteNotFound:
            warning = (
                f"Cannot add medium compound '{compound.id}' - metabolite not found in "
                f"extracellular compartment '{extracellular}'"
            )
            warnings.append(warning)
            logger.warning(warning)
        else:
            exchange_reaction = get_exchange_reaction(
                extracellular_metabolite, is_ec_model, consumption=True
            )

            # If someone already figured out the uptake rate for the compound, it's
            # likely more accurate than our assumptions, so keep it
            if exchange_reaction.id in model.medium:
                medium_mapping[exchange_reaction.id] = model.medium[
                    exchange_reaction.id
                ]
                continue

            if not extracellular_metabolite.formula:
                warning = (
                    f"No formula for metabolite '{extracellular_metabolite.id}', cannot"
                    f" check if it is a carbon source"
                )
                warnings.append(warning)
                logger.warning(warning)
                # If we don't know, it's most likely that the metabolite does not have a
                # higher uptake rate than a carbon source, so set the bound still to 10
                medium_mapping[exchange_reaction.id] = 10
            elif "C" in extracellular_metabolite.elements:
                # Limit the uptake rate for carbon sources to 10
                medium_mapping[exchange_reaction.id] = 10
            else:
                medium_mapping[exchange_reaction.id] = 1000

    # Apply the medium to the model, letting cobrapy deal with figuring out the correct
    # bounds to change
    model.medium = medium_mapping

    # Add all exchange reactions to operations, to make sure any changed bounds is
    # properly updated
    for reaction in model.exchanges:
        operations.append(
            {
                "operation": "modify",
                "type": "reaction",
                "id": reaction.id,
                "data": reaction_to_dict(reaction),
            }
        )

    return operations, warnings, errors


def apply_genotype(model, genotype_changes):
    """
    Apply genotype changes to a metabolic model.

    The genotype is first parsed by gnomic. Removed features (genes) are knocked
    out, while added features are added by looking up reaction definitions in
    ICE and adding those to the model.

    Parameters
    ----------
    model: cobra.Model
    genotype_changes: gnomic.Genotype
        A gnomic genotype object describing the strain modifications.

    Returns
    -------
    tuple (operations, warnings, errors)
        Operations is a list of model operations necessary to apply the medium
        to the model. See also the `Operations` schema.
        Warnings is a list of human-readable strings of potential issues.
        If errors is not an empty list, it was not possible to apply the
        genotype. Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    warnings = []
    errors = []

    # Apply feature operations
    for feature in genotype_changes.removed_features:
        feature_identifer = feature_id(feature)
        feature_lower = feature_identifer.lower()
        # Perform gene knockout. Use feature name as gene name
        try:
            # Some genotype descriptions wrongly use the protein names rather
            # than the gene names, for example, AdhE instead of adhE.
            # We want to be forgiving here and only compare lower case names.
            def compare_feature(gene):
                return (
                    gene.id == feature_identifer or gene.name.lower() == feature_lower
                )

            # We pick the first result. A fuzzy search on the name would be
            # useful in future.
            gene = model.genes.query(compare_feature)[0]
            gene.knock_out()
            operations.append({"operation": "knockout", "type": "gene", "id": gene.id})
        except IndexError:
            warning = (
                f"Cannot knockout gene '{feature_identifer}', not found in the model"
            )
            warnings.append(warning)
            logger.warning(warning)

    for feature in genotype_changes.added_features:
        feature_identifer = feature_id(feature)
        feature_lower = feature_identifer.lower()
        # Perform gene insertion unless the gene already exists in the model.

        def compare_feature(gene):
            return gene.id == feature_identifer or gene.name.lower() == feature_lower

        if model.genes.query(compare_feature):
            logger.info(
                f"Not adding gene '{feature_identifer}', "
                f"it already exists in the model."
            )
            continue

        try:
            heterologous = ice.get_reaction_equations(genotype=feature_identifer)
        except PartNotFound:
            warning = (
                f"Cannot add gene '{feature_identifer}', "
                f"no gene-protein-reaction rules were found on ICE."
            )
            warnings.append(warning)
            logger.warning(warning)
            continue

        for reaction_id, equation in heterologous.items():
            logger.info(
                f"Adding reaction '{reaction_id}' catalyzed by genetic part "
                f"'{feature_identifer}'."
            )
            if reaction_id in model.reactions:
                warning = (
                    f"Reaction {reaction_id} already exists in the "
                    f"model, removing and replacing it."
                )
                logger.warning(warning)
                warnings.append(warning)
                model.remove_reactions([reaction_id])
                operations.append(
                    {"operation": "remove", "type": "reaction", "id": reaction_id}
                )
            reaction = Reaction(reaction_id)
            reaction.gene_reaction_rule = feature_identifer
            model.add_reactions([reaction])

            # Before building the reaction's metabolites, keep track of the
            # existing ones to detect new metabolites added to the model.
            metabolites_before = set(model.metabolites)
            reaction.build_reaction_from_string(equation)
            new_metabolites = set(model.metabolites) - metabolites_before

            # Ensure all metabolites have a compartment. (Check all of the reaction's
            # metabolites, but presumably only new metabolites will not have a
            # compartment.)
            for metabolite in reaction.metabolites:
                if metabolite.compartment:
                    continue
                # Assume BiGG identifier convention and try to parse the compartment
                # id.
                try:
                    metabolite_id, compartment_id = parse_bigg_compartment(
                        metabolite.id, model
                    )
                except ValueError:
                    error = (
                        f"We cannot parse a compartment from heterologous metabolite "
                        f"'{metabolite.id}'."
                    )
                    errors.append(error)
                    logger.error(error)
                    continue
                except CompartmentNotFound as error:
                    error = (
                        f"Compartment {error.compartment_id} does not exist in the "
                        f"model, but that's what we understand metabolite "
                        f"{metabolite.id} to exist in."
                    )
                    errors.append(error)
                    logger.error(error)
                    continue
                logger.debug(
                    f"Setting compartment for metabolite {metabolite} to: "
                    f"{compartment_id}"
                )
                metabolite.compartment = compartment_id

            operations.append(
                {
                    "operation": "add",
                    "type": "reaction",
                    "data": reaction_to_dict(reaction),
                }
            )

            # We have to ensure that all new metabolites can leave the system. Create an
            # extracellular metabolite + exchange reaction (as opposed to just creating
            # a demand reaction for the intracellular metabolite) because if there are
            # metabolomics for this metabolite in a later step, our adapter logic always
            # assumes it exists in the 'e' compartment and that there exists an exchange
            # reaction.
            for metabolite in new_metabolites:
                # Create an extracellular version of the same metabolite
                if metabolite.compartment == "e":
                    logger.info(
                        f"Metabolite {metabolite} is already in the extracellular "
                        f"compartment; not creating transport/exchange reactions for it"
                    )
                    continue
                metabolite_id, compartment_id = parse_bigg_compartment(
                    metabolite.id, model
                )
                metabolite_e = Metabolite(
                    f"{metabolite_id}_e",
                    name=metabolite.name,
                    formula=metabolite.formula,
                    compartment="e",
                )

                # Create a transport reaction between the compartments
                transport_reaction = Reaction(
                    id=f"{metabolite_id.upper()}t", name=f"{metabolite.name} transport"
                )
                transport_reaction.bounds = Configuration().bounds
                transport_reaction.add_metabolites({metabolite: -1, metabolite_e: 1})
                model.add_reactions([transport_reaction])
                operations.append(
                    {
                        "operation": "add",
                        "type": "reaction",
                        "data": reaction_to_dict(transport_reaction),
                    }
                )

                # Create an exchange reaction for the extracellular metabolite so that
                # it may leave the system
                exchange_reaction = model.add_boundary(
                    metabolite_e, type="exchange", lb=0, ub=1000
                )
                operations.append(
                    {
                        "operation": "add",
                        "type": "reaction",
                        "data": reaction_to_dict(exchange_reaction),
                    }
                )

    return operations, warnings, errors


def apply_measurements(
    model,
    biomass_reaction,
    is_ec_model,
    fluxomics,
    metabolomics,
    proteomics,
    uptake_secretion_rates,
    molar_yields,
    growth_rate,
):
    """
    Apply omics measurements to a metabolic model.

    For each measured flux (production-rate / uptake-rate), constrain the model
    by forcing their upper and lower bounds to the measured values.

    Parameters
    ----------
    model: cobra.Model
    biomass_reaction: str
        The id of the biomass reaction in the given model.
    is_ec_model: bool
        A boolean indicating if the model is enzyme-constrained.
    fluxomics: list(dict)
        List of measurements matching the `Fluxomics` schema.
    metabolomics: list(dict)
        List of measurements matching the `Metabolomics` schema.
    proteomics: list(dict)
        List of measurements matching the `Proteomics` schema.
    uptake_secretion_rates: list(dict)
        List of measurements matching the `UptakeSecretionRates` schema.
    molar_yields: list(dict)
        List of measurements matching the `MolarYields` schema.
    growth_rate: dict
        Growth rate, matching the `GrowthRate` schema.

    Returns
    -------
    tuple (operations, warnings, errors)
        Operations is a list of model operations necessary to apply the
        measurements to the model. See also the `Operations` schema.
        Warnings is a list of human-readable strings of potential issues.
        If errors is not an empty list, it was not possible to apply the
        measurements.
        Errors then contains a list of string messages describing the
        problem(s).
    """
    operations = []
    warnings = []
    errors = []

    def bounds(measurement, uncertainty):
        """Return resolved bounds based on measurement and uncertainty"""
        if uncertainty:
            return (measurement - uncertainty, measurement + uncertainty)
        else:
            return (measurement, measurement)

    # First, improve the fluxomics dataset by minimizing the distance to a feasible
    # problem. If there is no objective constraint, skip minimization as it can yield
    # unreliable results.
    if growth_rate:
        growth_rate, fluxomics = minimize_distance(
            model, biomass_reaction, growth_rate, fluxomics
        )

    # If an enzyme constrained model with proteomics was supplied, flexibilize the
    # proteomics data and redefine the growth rate based on simulations.
    if growth_rate and proteomics and is_ec_model:
        growth_rate, proteomics, prot_warnings = flexibilize_proteomics(
            model, biomass_reaction, growth_rate, proteomics, uptake_secretion_rates
        )
        for warning in prot_warnings:
            warnings.append(warning)
            logger.warning(warning)

    # Constrain the model with the observed growth rate
    if growth_rate:
        reaction = model.reactions.get_by_id(biomass_reaction)
        reaction.bounds = bounds(growth_rate["measurement"], growth_rate["uncertainty"])
        operations.append(
            {
                "operation": "modify",
                "type": "reaction",
                "id": reaction.id,
                "data": reaction_to_dict(reaction),
            }
        )

    for measure in fluxomics:
        try:
            reaction = model.reactions.get_by_id(measure["identifier"])
        except KeyError:
            errors.append(
                f"Cannot find reaction '{measure['identifier']}' in the model"
            )
        else:
            reaction.bounds = bounds(measure["measurement"], measure["uncertainty"])
            operations.append(
                {
                    "operation": "modify",
                    "type": "reaction",
                    "id": reaction.id,
                    "data": reaction_to_dict(reaction),
                }
            )

    for metabolite in metabolomics:
        warning = (
            f"Cannot apply metabolomics measure for '{metabolite['identifier']}'; "
            f"feature has not yet been implemented"
        )
        warnings.append(warning)
        logger.warning(warning)

    for measure in proteomics:
        if is_ec_model:
            try:
                reaction = model.reactions.get_by_id(
                    f"prot_{measure['identifier']}_exchange"
                )
            except KeyError:
                warning = f"Cannot find protein '{measure['identifier']}' in the model"
                warnings.append(warning)
                logger.warning(warning)
            else:
                # measurement only modifies the upper bound (enzymes can be unsaturated)
                lb, ub = bounds(measure["measurement"], measure["uncertainty"])
                reaction.bounds = 0, ub
                operations.append(
                    {
                        "operation": "modify",
                        "type": "reaction",
                        "id": reaction.id,
                        "data": reaction_to_dict(reaction),
                    }
                )
        else:
            warning = (
                f"Cannot apply proteomics measurements for "
                f"non enzyme-constrained model {model.id}"
            )
            warnings.append(warning)
            logger.warning(warning)
            break

    for rate in uptake_secretion_rates:
        try:
            metabolite = find_metabolite(
                model, rate["identifier"], rate["namespace"], "e"
            )
        except MetaboliteNotFound as error:
            errors.append(str(error))
        else:
            exchange_reaction = get_exchange_reaction(
                metabolite, is_ec_model, consumption=rate["measurement"] < 0
            )
            lower_bound, upper_bound = bounds(rate["measurement"], rate["uncertainty"])

            # data is adjusted assuming a forward exchange reaction, i.e. x -->
            # (sign = -1), so if we instead actually have --> x, then multiply with -1
            # and flip lower bound and upper bound, to properly adjust for uncertainty,
            # e.g. if measurement = 3 and uncertainty = 0.3, then:
            # lb, ub = -1*(3 + 0.3), -1*(3 - 0.3) = -3.3, -2.7
            direction = exchange_reaction.metabolites[metabolite]
            if direction > 0:
                lower_bound, upper_bound = -1 * upper_bound, -1 * lower_bound
            exchange_reaction.bounds = lower_bound, upper_bound
            operations.append(
                {
                    "operation": "modify",
                    "type": "reaction",
                    "id": exchange_reaction.id,
                    "data": reaction_to_dict(exchange_reaction),
                }
            )

    for molar_yield in molar_yields:
        warning = (
            f"Cannot apply molar yield measurement for '"
            f"{molar_yield['product_identifier']}/{molar_yield['substrate_identifier']}"
            f"'; feature has not yet been implemented"
        )
        warnings.append(warning)
        logger.warning(warning)
    return operations, warnings, errors
