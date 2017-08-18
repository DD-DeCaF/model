import os
import json
import requests
from tqdm import tqdm

from cobra.io import read_sbml_model, write_sbml_model
from cameo import load_model

from model.adapter import add_prefix
from model.app import MODELS, MODEL_NAMESPACE
from model.settings import ID_MAPPER_API


LOCAL_MODELS = ['ecYeast7']
MODEL_METABOLITE_NAMESPACE = {
    'iJO1366': 'bigg.metabolite',
    'iMM904': 'bigg.metabolite',
    'iMM1415': 'bigg.metabolite',
    'iNJ661': 'bigg.metabolite',
    'iJN746': 'bigg.metabolite',
    'e_coli_core': 'bigg.metabolite',
    'ecYeast7': 'yeast7'
}

strip_compartment = {'yeast7': lambda mid: mid[:-2],
                     'bigg': lambda mid: mid[:-2]}


def sync_query_identifiers(object_ids, db_from, db_to):
    query = json.dumps({'ids': object_ids, 'dbFrom': db_from.lower(), 'dbTo': db_to.lower(), 'type': 'Metabolite'})
    r = requests.post(ID_MAPPER_API, data=query)
    return r.json()['ids']


def update_local_models(model_ids, model_store=None):
    """Update locally stored models.

    Annotate model metabolites with CHEBI identifiers and store them locally for easy access.
    :param model_ids: list, model identifiers
    :param model_store: path to directory where to store the processed models.
    """
    model_store = model_store or os.path.join(os.path.dirname(__file__), 'data/')
    for model_id in tqdm(model_ids):
        if model_id in LOCAL_MODELS:
            sbml_file = os.path.join(model_store, 'original', model_id + '.sbml.gz')
            model = read_sbml_model(sbml_file)
        else:
            model = load_model(model_id)
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
                if compound_id in {'s_0563', 'glc__D'} and db_name == 'CHEBI':
                    metabolite.annotation[db_name].append('CHEBI:42758')
                if metabolite_namespace not in metabolite.annotation:
                    metabolite.annotation[metabolite_namespace] = []
                metabolite.annotation[metabolite_namespace].append(compound_id)
        annotated_sbml_file = os.path.join(model_store, model_id + '.sbml.gz')

        write_sbml_model(model, annotated_sbml_file)

if '__main__' in __name__:
    for m_id in MODELS:
        update_local_models(m_id)
