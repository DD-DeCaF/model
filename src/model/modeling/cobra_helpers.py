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

from model.exceptions import MetaboliteNotFound, ReactionNotFound


logger = logging.getLogger(__name__)


def find_reaction(model, id, namespace):
    """
    Search for a given reaction in the model, also searching in annotations.

    Parameters
    ----------
    model: cobra.Model
    id: str
        The identifier of the reaction to find.
    namespace: str
        The miriam namespace identifier in which the given metabolite is
        registered. See https://www.ebi.ac.uk/miriam/main/collections

    Returns
    -------
    cobra.Reaction
        Returns the reaction object.

    Raises
    ------
    IndexError
        If multiple reaction are found for the given search query.
    ReactionNotFound
        If no reactions are found for the given parameters.
    """
    def query_fun(reaction):
        # Match namespace and identifiers case-insensitively
        for model_namespace, model_ids in reaction.annotation.items():
            if namespace.lower() != model_namespace.lower():
                return False
            if isinstance(model_ids, list):
                return id.lower() in [i.lower() for i in model_ids]
            else:
                return id.lower() == model_ids.lower()

    reactions = model.reactions.query(query_fun)
    if len(reactions) == 0:
        # No annotated result, try the default namespace (without confirming the
        # namespace id).
        try:
            return model.reactions.get_by_id(id)
        except KeyError:
            pass

        raise ReactionNotFound(
            f"Could not find reaction {id} in namespace {namespace} for "
            f"model {model.id}"
        )
    elif len(reactions) > 1:
        raise IndexError(f"Expected single reaction, found {reactions}")
    else:
        return reactions[0]


def find_metabolite(model, id, namespace, compartment):
    """
    Search for a given metabolite in the model, also searching in annotations.

    Parameters
    ----------
    model: cobra.Model
    id: str
        The identifier of the metabolite to find, e.g. "CHEBI:12965".
    namespace: str
        The miriam namespace identifier in which the given metabolite is
        registered. See https://www.ebi.ac.uk/miriam/main/collections
    compartment: str
        The compartment in which to look for the metabolite.

    Returns
    -------
    cobra.Metabolite
        Returns the metabolite object.

    Raises
    ------
    IndexError
        If multiple metabolites are found for the given search query.
    MetaboliteNotFound
        If no metabolites are found for the given parameters.
    """
    def query_fun(metabolite):
        if metabolite.compartment != compartment:
            return False

        # Try to find a case insensitive match for the namespace key
        for met_namespace in metabolite.annotation:
            if namespace.lower() == met_namespace.lower():
                identifier = metabolite.annotation[met_namespace]
                # Compare the identifier case insensitively as well
                # Annotations may contain a single id or a list of ids
                if isinstance(identifier, list):
                    return id.lower() in [i.lower() for i in identifier]
                else:
                    return id.lower() == identifier.lower()
        return False

    metabolites = model.metabolites.query(query_fun)
    if len(metabolites) == 0:
        # No annotated result, try the default namespace (without confirming the
        # namespace id).
        try:
            metabolite = model.metabolites.get_by_id(id)
            if metabolite.compartment == compartment:
                return metabolite
        except KeyError:
            pass

        raise MetaboliteNotFound(
            f"Could not find metabolite {id} in namespace {namespace} and "
            f"compartment {compartment} for model {model.id}"
        )
    elif len(metabolites) > 1:
        raise IndexError(f"Expected single metabolite, found {metabolites}")
    else:
        return metabolites[0]
