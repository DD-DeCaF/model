import logging
import pytest

from model.response import respond
from model.storage import Models

logging.disable(logging.CRITICAL)


@pytest.mark.asyncio
async def test_respond():
    message = {'to-return': ['fluxes', 'tmy', 'model', 'growth-rate', 'removed-reactions'], 'objectives': ['bigg:akg']}
    assert set((await respond(Models.get('iJO1366'), message)).keys()) == set(message['to-return'])
