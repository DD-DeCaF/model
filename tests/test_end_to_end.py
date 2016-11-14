import aiohttp
import pytest


@pytest.mark.asyncio
async def test_websocket():
    session = aiohttp.ClientSession()
    message_fluxes = {'to-return': ['fluxes']}
    message_modify = {
        'to-return': ['tmy', 'fluxes'],
        'objectives': ['chebi:17790'],
        'genotype-changes': ['+Aac'],
        'medium': [{'id': 'chebi:44080', 'concentration': 0.01}],
        'measurements': [{'id': 'chebi:44080', 'measurement': -15}],
    }
    async with session.ws_connect('http://localhost:8000/models/ECO') as ws:
        ws.send_json(message_fluxes)
        async for msg in ws:
            print(msg)
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
                    ws.send_json(message_modify)
            elif msg.type == aiohttp.WSMsgType.CLOSED:
                await ws.close()
                break
            elif msg.type == aiohttp.WSMsgType.ERROR:
                raise ws.exception()
