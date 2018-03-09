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
import os
import re

import requests
from cameo import load_model
from cobra.io import read_sbml_model, write_sbml_model
from tqdm import tqdm

from model.adapter import add_prefix
from model.constants import MODEL_NAMESPACE, MODELS
from model.settings import ID_MAPPER_API


LOCAL_MODELS = ['ecYeast7', 'ecYeast7_proteomics']
MODEL_METABOLITE_NAMESPACE = {
    'iJO1366': 'bigg.metabolite',
    'iMM904': 'bigg.metabolite',
    'iMM1415': 'bigg.metabolite',
    'iNJ661': 'bigg.metabolite',
    'iJN746': 'bigg.metabolite',
    'e_coli_core': 'bigg.metabolite',
    'ecYeast7': 'yeast7',
    'ecYeast7_proteomics': 'yeast7',
}

strip_compartment = {'yeast7': lambda mid: mid[:-2],
                     'bigg': lambda mid: mid[:-2]}

GLUCOSE = {'s_0563', 's_0565', 'glc__D', 's_0563', 's_0565', 's_0566', 's_0567', 's_0568', 's_1543'}


def sync_query_identifiers(object_ids, db_from, db_to):
    query = json.dumps({'ids': object_ids, 'dbFrom': db_from.lower(), 'dbTo': db_to.lower(), 'type': 'Metabolite'})
    r = requests.post(ID_MAPPER_API, data=query)
    return r.json()['ids']


def update_local_models(model_id, model_store=None):
    """Update locally stored models.

    Annotate model metabolites with CHEBI identifiers and store them locally for easy access.
    :param model_id: string, model identifier
    :param model_store: path to directory where to store the processed models.
    """
    model_store = model_store or '/io/model/data'
    if model_id in LOCAL_MODELS:
        sbml_file = os.path.join(model_store, 'original', model_id + '.sbml.gz')
        model = read_sbml_model(sbml_file)
    else:
        model = load_model(model_id)

    # annotate metabolites
    namespace = MODEL_NAMESPACE[model_id]
    metabolite_namespace = MODEL_METABOLITE_NAMESPACE[model_id]
    db_name = 'CHEBI'
    metabolites_missing_annotation = [m.id for m in model.metabolites if len(m.annotation.get(db_name, [])) < 1]
    model_xref = sync_query_identifiers(
        [strip_compartment[namespace](mid) for mid in metabolites_missing_annotation], namespace, db_name)
    for metabolite_id in metabolites_missing_annotation:
        compound_id = strip_compartment[namespace](metabolite_id)
        if compound_id in model_xref:
            metabolite = model.metabolites.get_by_id(metabolite_id)
            if db_name not in metabolite.annotation:
                metabolite.annotation[db_name] = []
            metabolite.annotation[db_name].extend(add_prefix(model_xref[compound_id], db_name))
            # TODO: For some reason, id-mapper doesn't make this link, add manually for now
            if compound_id in GLUCOSE and db_name == 'CHEBI':
                metabolite.annotation[db_name].append('CHEBI:42758')
            if metabolite_namespace not in metabolite.annotation:
                metabolite.annotation[metabolite_namespace] = []
            metabolite.annotation[metabolite_namespace].append(compound_id)

    # gecko protein exchanges
    db_name = 'uniprot'
    protein_exchanges = model.reactions.query(lambda rxn: re.match(r'^prot_.*_exchange$', rxn.id))
    for rxn in protein_exchanges:
        rxn.annotation[db_name] = [re.findall('^prot_(.*)_exchange$', rxn.id)[0]]
    write_sbml_model(model, os.path.join(model_store, model_id + '.sbml.gz'))

if '__main__' in __name__:
    for m_id in tqdm(MODELS):
        update_local_models(m_id)
