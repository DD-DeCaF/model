import aiohttp
import pytest
import requests

MESSAGE_FLUXES = {'to-return': ['fluxes'], 'measurements': [
    {'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0]},
    {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9]}]}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy'], 'objectives': ['chebi:17790'], 'request-id': 'requestid'}
MESSAGE_MODIFY = {
    'simulation-method': 'pfba',
    'to-return': ['tmy', 'fluxes'],
    'objectives': ['chebi:17790'],
    'genotype-changes': ['+Aac'],
    'measurements': [{'unit': 'mmol', 'name': 'aldehydo-D-glucose', 'id': 'chebi:42758', 'measurements': [-9.0]},
                     {'unit': 'mmol', 'name': 'ethanol', 'id': 'chebi:16236', 'measurements': [5.0, 4.8, 5.2, 4.9]}],
}

URL = 'http://localhost:8000/models/'
MODEL_OPTIONS_URL = 'http://localhost:8000/model-options/ECOLX'
URL_MAPS = 'http://localhost:8000/maps'
WS_URL = 'http://localhost:8000/wsmodels/'


def test_http():
    r = requests.get(MODEL_OPTIONS_URL)
    r.raise_for_status()
    r = requests.get(URL_MAPS)
    r.raise_for_status()
    for message in [MESSAGE_MODIFY, MESSAGE_FLUXES]:
        r = requests.post(URL + 'iJO1366', json={'message': message})
        r.raise_for_status()
        assert set(r.json().keys()) == set(message['to-return'] + ['model-id'])
        assert r.json()['fluxes']['EX_glc__D_e'] == -9.0
        assert abs(r.json()['fluxes']['EX_etoh_e'] - 4.64) < 0.001  # lower bound
    r = requests.post(URL + 'wrong_id', json={'message': {}})
    assert r.status_code == 404
    r = requests.post(URL + 'iJO1366', json={})
    assert r.status_code == 400


@pytest.mark.asyncio
async def test_websocket():
    r = requests.post(URL + 'iJO1366', json={'message': MESSAGE_MODIFY})
    r.raise_for_status()
    model_id = r.json()['model-id']
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
