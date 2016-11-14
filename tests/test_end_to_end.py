import aiohttp
import pytest
import requests

MESSAGE_FLUXES = {'to-return': ['fluxes']}
MESSAGE_MODIFY = {
    'to-return': ['tmy', 'fluxes'],
    'objectives': ['chebi:17790'],
    'genotype-changes': ['+Aac'],
    'medium': [{'id': 'chebi:44080', 'concentration': 0.01}],
    'measurements': [{'id': 'chebi:44080', 'measurement': -15}],
}


def test_http():
    URL = 'http://localhost:8000/models/'
    for message in [MESSAGE_MODIFY, MESSAGE_FLUXES]:
        r = requests.post(URL + 'iJO1366', json={'message': message})
        r.raise_for_status()
        assert set(r.json().keys()) == set(message['to-return'])
    r = requests.post(URL + 'wrong_id', json={'message': {}})
    assert r.status_code == 404
    r = requests.post(URL + 'iJO1366', json={})
    assert r.status_code == 400


@pytest.mark.asyncio
async def test_websocket():
    session = aiohttp.ClientSession()
    async with session.ws_connect('http://localhost:8000/wsmodels/iJO1366') as ws:
        ws.send_json(MESSAGE_FLUXES)
        async for msg in ws:
            if msg.type == aiohttp.WSMsgType.TEXT:
                if msg.data == 'close':
                    await ws.close()
                    break
                if 'tmy' in msg.json():
                    assert msg.json()['fluxes']
                    assert msg.json()['tmy']
                    await ws.close()
                    break
                else:
                    ws.send_json(MESSAGE_MODIFY)
            elif msg.type == aiohttp.WSMsgType.CLOSED:
                await ws.close()
                break
            elif msg.type == aiohttp.WSMsgType.ERROR:
                raise ws.exception()
