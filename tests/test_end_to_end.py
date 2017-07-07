import aiohttp
import pytest
import requests

MEASUREMENTS = [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0],
                 'type': 'compound'},
                {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9],
                 'type': 'compound'},
                {'id': 'PFK', 'measurements': [5, 5], 'type': 'reaction'}]
MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': MEASUREMENTS}
MESSAGE_FLUXES_INFEASIBLE = {'to-return': ['fluxes'], 'measurements': [
    {'id': 'ATPM', 'measurements': [100, 100], 'type': 'reaction'}]}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy'], 'objectives': ['chebi:17790'], 'request-id': 'requestid'}
MESSAGE_MODIFY = {
    'simulation-method': 'pfba',
    "reactions-add": [
        "MNXR81835", "MNXR83321"
    ],
    'to-return': ['tmy', 'fluxes', 'model', 'added-reactions', 'removed-reactions'],
    'objectives': ['chebi:17790'],
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
            changes = response.json()['model']['notes']['changes']
            assert 'PTA2' in {rxn['id'] for rxn in changes['removed']['reactions']}
            assert 'EX_etoh_e' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'PFK' in {rxn['id'] for rxn in changes['measured']['reactions']}
            assert 'b2297' in {rxn['id'] for rxn in changes['removed']['genes']}
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
    session = aiohttp.ClientSession()
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
