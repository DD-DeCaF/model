# Copyright 2018 Novo Nordisk Foundation Center for Biosustainability, DTU.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging

import pytest

from model.response import respond
from model.storage import Models


logging.disable(logging.CRITICAL)


@pytest.mark.asyncio
async def test_respond():
    message = {'to-return': ['fluxes', 'tmy', 'model', 'growth-rate', 'removed-reactions'], 'objectives': ['bigg:akg']}
    assert set((await respond(Models.get('iJO1366'), message)).keys()) == set(message['to-return'])
