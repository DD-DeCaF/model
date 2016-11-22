import pytest
import random
import json
from cobra.io.json import _to_dict
from model.app import existing_metabolite, NoIDMapping, restore_model, find_in_memory, product_reaction_variable, \
    phase_plane_to_dict, new_features_identifiers, apply_reactions_knockouts, respond, save_changes_to_db, \
    key_from_model_info, GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, convert_mg_to_mmol, convert_measurements_to_mmol, \
    modify_model, restore_from_db
from driven.generic.adapter import full_genotype


def almost_equal(a, b):
    return abs(a - b) < 10e-6


def test_mg_to_mmol():
    assert almost_equal(convert_mg_to_mmol(34.5, 58.4), 0.59075342)
    assert almost_equal(convert_mg_to_mmol(18, 18), 1)


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
    message = {
        GENOTYPE_CHANGES: ['+Aac'],
        MEASUREMENTS: [{'id': 'chebi:44080', 'concentration': 0.01}],
        MEDIUM: [{'id': 'chebi:44080', 'measurement': -15, 'unit': 'mg', 'name': 'glucose'}],
    }
    model = await modify_model(message, (await restore_model(model_id)).copy())
    db_key = await save_changes_to_db(model, model_id, message)
    assert len(_to_dict(await restore_from_db(db_key))['reactions']) == len(_to_dict(model)['reactions'])


def test_existing_metabolite():
    ecoli = find_in_memory('iJO1366')
    assert existing_metabolite(ecoli, 'chebi:17790') == existing_metabolite(ecoli, 'bigg:meoh')
    assert existing_metabolite(ecoli, 'bigg:succ').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        existing_metabolite(ecoli, 'wrong_id')


def test_product_reaction_variable():
    ecoli = find_in_memory('iJO1366')
    assert product_reaction_variable(ecoli, 'bigg:akg').id == 'EX_akg_e'
    assert product_reaction_variable(ecoli, 'bigg:e4p') is None


def test_phase_plane_to_dict():
    ecoli = find_in_memory('iJO1366')
    result = phase_plane_to_dict(ecoli, 'bigg:glu__L')
    assert set(result.keys()) == {'EX_glu__L_e', 'objective_lower_bound', 'objective_upper_bound'}
    assert len(set([len(v) for v in result.values()])) == 1
    assert phase_plane_to_dict(ecoli, 'bigg:g3p') == {}


def test_new_features_identifiers():
    changes = full_genotype(['-A -B +promoter.C:#D:#E:terminator.F', '+G', '+Y -H'])
    result = new_features_identifiers(changes)
    assert set(result) == {'C', 'D', 'E', 'F', 'G', 'Y'}


def test_respond():
    message = {'to-return': ['fluxes', 'tmy', 'model'], 'objectives': ['bigg:akg']}
    assert set(respond(message, find_in_memory('iJO1366')).keys()) == set(message['to-return'])


@pytest.mark.asyncio
async def test_apply_reactions_knockouts():
    ecoli = find_in_memory('iJO1366')
    result = (await apply_reactions_knockouts(ecoli, ['GLUDy', '3HAD160', 'GLUDy'])).model
    assert result.reactions.GLUDy.lower_bound == result.reactions.GLUDy.upper_bound == 0


def test_convert_measurements_to_mmol():
    ecoli = find_in_memory('iJO1366')
    measurements = [{'id': 'chebi:17790', 'measurement': 32.04186, 'unit': 'mg'}]
    assert convert_measurements_to_mmol(measurements, ecoli) == [{'id': 'chebi:17790', 'measurement': 1.0, 'unit': 'mmol'}]
