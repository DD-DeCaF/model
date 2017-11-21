import logging
import random
import pytest

from cobra.io import model_to_dict

from model.constants import GENOTYPE_CHANGES, MEASUREMENTS, MEDIUM
from model.operations import modify_model
from model.storage import (key_from_model_info, save_changes_to_db, find_in_memory, restore_model,
                           restore_from_db)

logging.disable(logging.CRITICAL)

def test_key_from_model_info():
    model_id = 'model_id'
    for _ in range(50):
        values1 = random.sample(range(20), 10)
        values2 = random.sample(range(20), 10)
        values3 = random.shuffle(values2[:])
        message1 = {MEASUREMENTS: values1, GENOTYPE_CHANGES: values2}
        message2 = {MEASUREMENTS: values2}
        message3 = {MEASUREMENTS: values3, GENOTYPE_CHANGES: values2}
        message4 = {GENOTYPE_CHANGES: values2, MEASUREMENTS: values1}
        assert key_from_model_info(model_id, message1) != key_from_model_info(model_id, message2)
        assert key_from_model_info(model_id, message1) != key_from_model_info(model_id, message3)
        assert key_from_model_info(model_id, message1) == key_from_model_info(model_id, message4)


@pytest.mark.asyncio
async def test_save_and_restore():
    model_id = 'e_coli_core'
    await save_changes_to_db(find_in_memory(model_id), model_id, {})
    message = {
        GENOTYPE_CHANGES: [
            '-aceA -sucCD -pykA -pykF -pta +promoter.BBa_J23100:#AB326105:#NP_600058:terminator.BBa_0010'],
        MEASUREMENTS: [{'id': 'chebi:44080', 'measurements': [0.01, float('nan')], 'type': 'compound'}],
        MEDIUM: [{'concentration': 27.0, 'id': 'chebi:42758'}, {'concentration': 6.0, 'id': 'chebi:16015'},
                 {'concentration': 1.6, 'id': 'chebi:30808'}, {'concentration': 2.0, 'id': 'chebi:35696'},
                 {'concentration': 1.0, 'id': 'chebi:49553'}, {'concentration': 2.0, 'id': 'chebi:49976'},
                 {'concentration': 0.5, 'id': 'chebi:33118'}, {'concentration': 3.5, 'id': 'chebi:63036'},
                 {'concentration': 5.0, 'id': 'chebi:131527'}, {'concentration': 3.5, 'id': 'chebi:63051'},
                 {'concentration': 0.5, 'id': 'chebi:31795'}, {'concentration': 0.005, 'id': 'chebi:3312'},
                 {'concentration': 0.005, 'id': 'chebi:49105'}, {'concentration': 300.0, 'id': 'chebi:42758'},
                 {'concentration': 9.0, 'id': 'chebi:16015'}, {'concentration': 52.5, 'id': 'chebi:62946'}],
    }
    model = await modify_model(message, (await restore_model(model_id)).copy())
    db_key = await save_changes_to_db(model, model_id, message)
    restored = model_to_dict(await restore_from_db(db_key))
    original = model_to_dict(model)
    assert set([i['id'] for i in restored['reactions']]) == set([i['id'] for i in original['reactions']])
    assert set([i['id'] for i in restored['genes']]) == set([i['id'] for i in original['genes']])
    assert set([i['id'] for i in restored['metabolites']]) == set([i['id'] for i in original['metabolites']])


@pytest.mark.asyncio
async def test_model_immutability():
    """Changes on restored models must not affect cache"""
    model = (await restore_model('e_coli_core')).copy()
    model.notes['test'] = 'test'
    restored_model = (await restore_model('e_coli_core')).copy()
    restored_model.notes['test'] = 'different'
    assert model.notes['test'] == 'test'
