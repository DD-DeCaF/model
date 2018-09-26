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

from cobra.io import read_sbml_model

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


@pytest.fixture(scope="session")
def models():
    """
    Provide loaded and instantiated cobrapy models for test usage. This fixture ensures models are loaded only once per
    test session, but please use the function-scoped fixtures below to be able to make revertable modifications to the
    models."""
    return {
        'e_coli_core': read_sbml_model('tests/data/e_coli_core.sbml.gz'),
        'iJO1366': read_sbml_model('tests/data/iJO1366.sbml.gz'),
    }


@pytest.fixture(scope="function")
def e_coli_core(models):
    """
    Provide the e_coli_core model in a context manager, so that modifications are not persisted beyond the scope of the
    test function. This model is fairly small and should be preferred in test cases where possible.
    """
    with models['e_coli_core'] as model:
        yield model


@pytest.fixture(scope="function")
def iJO1366(models):
    """
    Provide the iJO1366 model in a context manager, so that modifications are not persisted beyond the scope of the
    test function.
    """
    with models['iJO1366'] as model:
        yield model
