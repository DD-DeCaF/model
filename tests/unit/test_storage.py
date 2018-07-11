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

import random

from cobra.io.dict import model_to_dict

from model import storage
from model.constants import GENOTYPE_CHANGES, MEASUREMENTS, MEDIUM
from model.ice_client import ICE
from model.operations import modify_model


def test_key_from_model_info(app):
    model_id = 'model_id'
    for _ in range(50):
        values1 = random.sample(range(20), 10)
        values2 = random.sample(range(20), 10)
        values3 = random.shuffle(values2[:])
        message1 = {MEASUREMENTS: values1, GENOTYPE_CHANGES: values2}
        message2 = {MEASUREMENTS: values2}
        message3 = {MEASUREMENTS: values3, GENOTYPE_CHANGES: values2}
        message4 = {GENOTYPE_CHANGES: values2, MEASUREMENTS: values1}
        assert storage._changes_key(model_id, message1) != storage._changes_key(model_id, message2)
        assert storage._changes_key(model_id, message1) != storage._changes_key(model_id, message3)
        assert storage._changes_key(model_id, message1) == storage._changes_key(model_id, message4)


def test_save_and_restore(monkeypatch, app):
    # Disable GPR queries for efficiency
    monkeypatch.setattr(ICE, 'get_reaction_equations', lambda self, genotype: {})

    model_id = 'e_coli_core'
    storage.save_changes(storage.get(model_id).model, {})
    message = {
        GENOTYPE_CHANGES: [
            '-aceA -sucCD -pykA -pykF -pta +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'],
        MEASUREMENTS: [{'id': 'chebi:44080', 'measurements': [0.01, float('nan')], 'type': 'compound'},
                       {'measurements': [0.01, float('nan')], 'type': 'growth-rate'}],
        MEDIUM: [{'concentration': 27.0, 'id': 'chebi:42758'}, {'concentration': 6.0, 'id': 'chebi:16015'},
                 {'concentration': 1.6, 'id': 'chebi:30808'}, {'concentration': 2.0, 'id': 'chebi:35696'},
                 {'concentration': 1.0, 'id': 'chebi:49553'}, {'concentration': 2.0, 'id': 'chebi:49976'},
                 {'concentration': 0.5, 'id': 'chebi:33118'}, {'concentration': 3.5, 'id': 'chebi:63036'},
                 {'concentration': 5.0, 'id': 'chebi:131527'}, {'concentration': 3.5, 'id': 'chebi:63051'},
                 {'concentration': 0.5, 'id': 'chebi:31795'}, {'concentration': 0.005, 'id': 'chebi:3312'},
                 {'concentration': 0.005, 'id': 'chebi:49105'}, {'concentration': 300.0, 'id': 'chebi:42758'},
                 {'concentration': 9.0, 'id': 'chebi:16015'}, {'concentration': 52.5, 'id': 'chebi:62946'}],
    }
    model = modify_model(message, storage.get(model_id).model.copy())
    db_key = storage.save_changes(model, message)
    restored = model_to_dict(storage.restore_from_key(db_key))
    original = model_to_dict(model)
    assert set([i['id'] for i in restored['reactions']]) == set([i['id'] for i in original['reactions']])
    assert set([i['id'] for i in restored['genes']]) == set([i['id'] for i in original['genes']])
    assert set([i['id'] for i in restored['metabolites']]) == set([i['id'] for i in original['metabolites']])


def test_model_immutability(app):
    """Changes on restored models must not affect cache"""
    model = storage.get('e_coli_core').model.copy()
    model.notes['test'] = 'test'
    restored_model = storage.get('e_coli_core').model
    restored_model.notes['test'] = 'different'
    assert model.notes['test'] == 'test'
