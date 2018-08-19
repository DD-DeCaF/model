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

import pytest

from model import storage
from model.app import app as app_
from model.app import init_app


@pytest.fixture(scope="session")
def app():
    """Provide the initialized Flask app"""
    init_app(app_, None)
    return app_


@pytest.fixture(scope="session")
def client(app):
    """Provide a Flask test client"""
    with app.test_client() as client:
        yield client


@pytest.fixture(scope="function")
def e_coli_core():
    """
    Provide a modifiable copy of the e_coli_core cobrapy model.  This model is fairly small and should be preferred in
    test cases where possible.
    """
    with storage.get('e_coli_core').model as model:
        yield model


@pytest.fixture(scope="function")
def iJO1366():
    """Provide a modifiable copy of the iJO1366 cobrapy model"""
    with storage.get('iJO1366').model as model:
        yield model


@pytest.fixture(scope="function")
def iMM904():
    """Provide a modifiable copy of the iMM904 cobrapy model"""
    with storage.get('iMM904').model as model:
        yield model
