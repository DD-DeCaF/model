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

from simulations.exceptions import (
    CompartmentNotFound,
    MetaboliteNotFound,
    ReactionNotFound,
)


logger = logging.getLogger(__name__)


def find_reaction(model, id, namespace):
    """
    Search a model for a given reaction, also looking in annotations.

    Parameters
    ----------
    model: cobra.Model
    id: str
        The identifier of the reaction to find. The comparison is made case
        insensitively.
    namespace: str
        The miriam namespace identifier in which the given metabolite is
        registered. See https://www.ebi.ac.uk/miriam/main/collections
        The comparison is made case insensitively.

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
        return _query_item(reaction, id, namespace)

    reactions = model.reactions.query(query_fun)
    if len(reactions) == 0:
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
    Search a model for a given metabolite, also looking in annotations.

    Parameters
    ----------
    model: cobra.Model
    id: str
        The identifier of the metabolite to find, e.g. "CHEBI:12965". The
        comparison is made case insensitively. If a match is not found, another
        attempt will be made with the compartment id appended, e.g., looking for
        both "o2" and "o2_e". However, the caller should ideally pass the full
        metabolite identifier.
    namespace: str
        The miriam namespace identifier in which the given metabolite is
        registered. See https://www.ebi.ac.uk/miriam/main/collections
        The comparison is made case insensitively.
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

        result = _query_item(metabolite, id, namespace)
        if result:
            return result

        # If the original query fails, retry with the compartment id appended
        # to the identifier (a regular convenation with BiGG metabolites, but
        # may also be the case in other namespaces).
        return _query_item(metabolite, f"{id}_{compartment}", namespace)

    metabolites = model.metabolites.query(query_fun)
    if len(metabolites) == 0:
        raise MetaboliteNotFound(
            f"Could not find metabolite {id} or {id}_{compartment} in "
            f"namespace {namespace} and compartment {compartment} for model "
            f"{model.id}"
        )
    elif len(metabolites) > 1:
        raise IndexError(f"Expected single metabolite, found {metabolites}")
    else:
        return metabolites[0]


def _query_item(item, query_id, query_namespace):
    """
    Check if the given cobra collection item matches the query arguments.

    Parameters
    ----------
    item: cobra.Reaction or cobra.Metabolite
    query_id: str
        The identifier to compare. The comparison is made case insensitively.
    query_namespace: str
        The miriam namespace identifier in which the given metabolite is
        registered. See https://www.ebi.ac.uk/miriam/main/collections
        The comparison is made case insensitively.

    Returns
    -------
    bool
        True if the given id exists in the default namespace, or in the model
        annotations by the queried namespace, otherwise False.
    """
    # Try the default identifiers (without confirming the namespace)
    if query_id.lower() == item.id.lower():
        return True

    # Otherwise, try to find a case insensitive match for the namespace key
    for namespace in item.annotation:
        if query_namespace.lower() == namespace.lower():
            annotation = item.annotation[namespace]
            # Compare the identifier case insensitively as well
            # Annotations may contain a single id or a list of ids
            if isinstance(annotation, list):
                if query_id.lower() in [i.lower() for i in annotation]:
                    return True
            else:
                if query_id.lower() == annotation.lower():
                    return True
    return False


def parse_bigg_compartment(metabolite_id, model):
    """
    Parse the compartment ID of the given metabolite identifier.

    Parameters
    ----------
    metabolite_id: str
        A metabolite identifier, expected to follow the BiGG convention of appending an
        underscore followed by the compartment id. For example: "h2o_c"
    model: cobra.Model
        The model the metabolite belongs to. Used to confirm that the parsed compartment
        id exists in the model.

    Returns
    -------
    tuple(str, str)
        First string: The metabolite identifier without its compartment postfix.
        Second string: The compartment identifier.

    Raises
    ------
    ValueError
        If the given identifier does not match the expected BiGG format, i.e., does not
        contain an underscore.
    CompartmentNotFound
        If the compartment identifier parsed out of the metabolite does not exist in the
        model.
    """
    try:
        metabolite_id, compartment_id = metabolite_id.rsplit("_", 1)
    except ValueError as error:
        raise ValueError(
            f"The identifier {metabolite_id} has no valid BiGG compartment suffix."
        ) from error
    if compartment_id not in model.compartments:
        raise CompartmentNotFound(
            f"Compartment {compartment_id} does not exist in model {model.id}",
            compartment_id,
        )
    return metabolite_id, compartment_id
