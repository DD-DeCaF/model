import aiohttp
import pytest
import requests

MESSAGE_FLUXES = {'to-return': ['fluxes']}
MESSAGE_TMY_FLUXES = {'to-return': ['fluxes', 'tmy'], 'objectives': ['chebi:17790']}
MESSAGE_MODIFY = {
    'to-return': ['tmy', 'fluxes'],
    'objectives': ['chebi:17790'],
    'genotype-changes': ['+Aac'],
}

URL = 'http://localhost:8000/models/'
WS_URL = 'http://localhost:8000/wsmodels/'


def test_http():
    for message in [MESSAGE_MODIFY, MESSAGE_FLUXES]:
        r = requests.post(URL + 'iJO1366', json={'message': message})
        r.raise_for_status()
        assert set(r.json().keys()) == set(message['to-return'] + ['model-id'])
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
