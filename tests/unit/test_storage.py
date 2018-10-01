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

import requests
from cobra import Model

from model import storage


class MockResponse:
    status_code = 200

    def json(self):
        return {
            'model_serialized': {
                'version': "1",
                'id': "foo",
                'genes': [],
                'reactions': [],
                'metabolites': [],
                'compartments': {},
            },
            'organism_id': "bar",
            'default_biomass_reaction': "baz",
        }

    def raise_for_status(self):
        pass


def test_get_model(monkeypatch):
    monkeypatch.setattr(requests, 'get', lambda url: MockResponse())
    assert type(storage.get(10).model) == Model
