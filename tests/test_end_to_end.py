import aiohttp
from copy import deepcopy
import pytest
import requests
import json
import jsonpatch
import unittest

MEASUREMENTS = [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0],
                 'type': 'compound'},
                {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9],
                 'type': 'compound'},
                {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'},
                {'id': 'BAD_ID', 'measurements': [5, 5], 'type': 'reaction', 'db_name': 'bigg.reaction'}]
MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': MEASUREMENTS}
MESSAGE_FLUXES_INFEASIBLE = {'to-return': ['fluxes'], 'measurements': [
    {'id': 'ATPM', 'measurements': [100, 100], 'type': 'reaction', 'db_name': 'bigg.reaction'}]}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy'], 'objectives': ['chebi:17790'], 'request-id': 'requestid'}
MESSAGE_MODIFY = {
    'simulation-method': 'pfba',
    'reactions-add': [
        {'id': 'MNXR69355', 'metabolites': None},
        {'id': 'MNXR83321', 'metabolites': None},
        {'id': 'MagicCarrot', 'metabolites': {'glc__D_c': -1, 'caro_c': 1}}
    ],
    'to-return': ['tmy', 'fluxes', 'model', 'added-reactions', 'removed-reactions'],
    'objectives': ['chebi:17790', 'chebi:17579'],
    'genotype-changes': ['+Aac', '-pta'],
    'measurements': MEASUREMENTS,
}

URL = 'http://localhost:8000/models/'
MODEL_OPTIONS_URL = 'http://localhost:8000/model-options/ECOLX'
URL_MAPS = 'http://localhost:8000/maps'
WS_URL = 'http://localhost:8000/wsmodels/'


def test_http():
    response = requests.get(MODEL_OPTIONS_URL)
    response.raise_for_status()
    response = requests.get(URL_MAPS)
    response.raise_for_status()
    for query, message in {'modify': MESSAGE_MODIFY, 'fluxes': MESSAGE_FLUXES}.items():
        response = requests.post(URL + 'iJO1366', json={'message': message})
        response.raise_for_status()
        assert set(response.json().keys()) == set(message['to-return'] + ['model-id'])
        assert response.json()['fluxes']['EX_glc__D_e'] == -9.0
        assert abs(response.json()['fluxes']['EX_etoh_e'] - 4.64) < 0.001  # lower bound
        assert response.json()['fluxes']['PFK'] == 5
        if query == 'modify':
            tmy = response.json()['tmy']
            changes = response.json()['model']['notes']['changes']
            assert sum(tmy['chebi:17579']['DM_caro_e']) > 10.
            assert 'PTA2' in {rxn['id'] for rxn in changes['removed']['reactions']}
            assert 'EX_etoh_e' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'PFK' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'b2297' in {rxn['id'] for rxn in changes['removed']['genes']}
            assert 'BAD_ID' in {rxn['id'] for rxn in changes['measured-missing']['reactions']}
    response = requests.post(URL + 'iJO1366', json={'message': MESSAGE_FLUXES_INFEASIBLE})
    assert response.json()['fluxes']['ATPM'] == 100
    response = requests.post(URL + 'wrong_id', json={'message': {}})
    assert response.status_code == 404
    response = requests.post(URL + 'iJO1366', json={})
    assert response.status_code == 400


@pytest.mark.asyncio
async def test_websocket():
    response = requests.post(URL + 'iJO1366', json={'message': MESSAGE_MODIFY})
    response.raise_for_status()
    model_id = response.json()['model-id']
    async with aiohttp.ClientSession() as session:
        async with session.ws_connect(WS_URL + model_id) as ws:
            ws.send_json(MESSAGE_TMY_FLUXES)
            async for msg in ws:
                if msg.type == aiohttp.WSMsgType.TEXT:
                    assert msg.json()['fluxes']
                    assert msg.json()['tmy']
                elif msg.type == aiohttp.WSMsgType.ERROR:
                    raise ws.exception()
                await ws.close()
                break


@pytest.mark.asyncio
async def test_diff_model_api():
    original_message = deepcopy(MESSAGE_MODIFY)
    original_message['to-return'] = ["fluxes", "model", "added-reactions", "removed-reactions"]
    original_response = requests.post('{url}{model}'.format(url=URL, model='iJO1366'), json={'message': original_message})
    original_response.raise_for_status()
    original_data = original_response.json()

    model_response = requests.get('{url}{model}'.format(url='http://localhost:8000/v1/models/', model='iJO1366'))
    model_response.raise_for_status()

    diff_response = requests.post('{url}{model}'.format(url='http://localhost:8000/v1/models/', model='iJO1366'), json={'message': original_message})
    patch = jsonpatch.JsonPatch.from_string(diff_response.text)
    result = patch.apply(model_response.json())

    del result['model']['notes']['changes']
    del original_data['model']['notes']['changes']
    tc = unittest.TestCase()
    tc.assertDictEqual(result['model'], original_data['model'])
