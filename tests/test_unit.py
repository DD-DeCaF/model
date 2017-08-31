import pytest
import random
from copy import deepcopy
from cobra.io import model_to_dict, read_sbml_model
from tempfile import mkdtemp
import os

from model.app import (NoIDMapping, restore_model, find_in_memory, product_reaction_variable,
                       phase_plane_to_dict, new_features_identifiers, apply_reactions_knockouts, respond,
                       save_changes_to_db,
                       key_from_model_info, GENOTYPE_CHANGES, MEDIUM, MEASUREMENTS, convert_mg_to_mmol,
                       convert_measurements_to_mmol,
                       modify_model, restore_from_db, add_reactions, EMPTY_CHANGES)
from model.adapter import full_genotype, get_unique_metabolite
from model.update_models import update_local_models


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
async def test_b_number():
    model_id = 'iJO1366'
    message = {
        GENOTYPE_CHANGES: ['-b3067,-b3172,-b1084'],
    }
    model = await modify_model(message, (await restore_model(model_id)).copy())
    assert not model.genes.b3172.functional


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


@pytest.mark.asyncio
async def test_existing_metabolite():
    ecoli = find_in_memory('iJO1366')
    assert get_unique_metabolite(ecoli, 'chebi:17790') == get_unique_metabolite(
        ecoli, 'meoh', db_name='bigg.metabolite')
    assert get_unique_metabolite(ecoli, 'succ', db_name='bigg.metabolite').formula == 'C4H4O4'
    with pytest.raises(NoIDMapping):
        await get_unique_metabolite(ecoli, 'wrong_id')


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


@pytest.mark.asyncio
async def test_respond():
    message = {'to-return': ['fluxes', 'tmy', 'model', 'growth-rate', 'removed-reactions'], 'objectives': ['bigg:akg']}
    assert set((await respond(message, find_in_memory('iJO1366'))).keys()) == set(message['to-return'])


@pytest.mark.asyncio
async def test_reactions_knockouts():
    ecoli_original = find_in_memory('iJO1366').copy()
    ecoli = ecoli_original.copy()
    ecoli.notes['changes'] = deepcopy(EMPTY_CHANGES)
    reaction_ids = {'GLUDy', 'GLUDy', '3HAD160'}
    GLUDy_upper_bound = ecoli.reactions.get_by_id('GLUDy').upper_bound
    assert GLUDy_upper_bound != 0
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == reaction_ids
    assert ecoli.reactions.get_by_id('GLUDy').upper_bound == 0
    reaction_ids = reaction_ids - {'GLUDy'}
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == {'3HAD160'}
    assert GLUDy_upper_bound == ecoli.reactions.get_by_id('GLUDy').upper_bound
    reaction_ids = reaction_ids - {'3HAD160'}
    ecoli = await apply_reactions_knockouts(ecoli, list(reaction_ids))
    assert set([i['id'] for i in ecoli.notes['changes']['removed']['reactions']]) == set()
    assert almost_equal(ecoli.optimize().objective_value, ecoli_original.optimize().objective_value)


@pytest.mark.asyncio
async def test_convert_measurements_to_mmol():
    ecoli = find_in_memory('iJO1366').copy()
    measurements = [{'id': 'chebi:17790', 'measurements': [32.04186], 'unit': 'mg', 'type': 'compound'}]
    assert convert_measurements_to_mmol(measurements, ecoli) == [
        {'id': 'chebi:17790', 'measurements': [1.0], 'unit': 'mmol', 'type': 'compound'}]


def test_add_reactions():
    ecoli = find_in_memory('iJO1366').copy()
    keys = ('id', 'name', 'metabolites', 'lower_bound', 'upper_bound', 'gene_reaction_rule')
    info = ('newid', '', {}, 0, 0, '')
    changes = [dict(zip(keys, info))] * 2
    add_reactions(ecoli, changes)


def test_update_models():
    tempdir = mkdtemp()
    update_local_models('e_coli_core', tempdir)
    model = read_sbml_model(os.path.join(tempdir, 'e_coli_core.sbml.gz'))
    glucose = get_unique_metabolite(model, 'CHEBI:42758')
    assert glucose.id == 'glc__D_e'
    assert glucose.annotation['bigg.metabolite'] == 'glc__D'
