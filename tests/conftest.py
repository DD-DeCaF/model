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


@pytest.fixture(scope="session")
def models():
    """
    Preloads the storage module with the two test models. This fixture ensures models are loaded locally, and only once
    per test session. Use the identifiers 'test_e_coli_core' and 'test_iJO1366' for endpoint tests.
    For unit tests, consider using the function-scoped fixtures below to be able to make revertable modifications to the
    models.
    """
    model = read_sbml_model('tests/data/e_coli_core.sbml.gz')
    storage._MODELS['test_e_coli_core'] = storage.ModelWrapper(model, None, "Escherichia coli",
                                                               'BIOMASS_Ecoli_core_w_GAM')
    model = read_sbml_model('tests/data/e_coli_core.sbml.gz')
    storage._MODELS['test_e_coli_core_proprietary'] = storage.ModelWrapper(model, 1, "Escherichia coli",
                                                               'BIOMASS_Ecoli_core_w_GAM')
    model = read_sbml_model('tests/data/iJO1366.sbml.gz')
    storage._MODELS['test_iJO1366'] = storage.ModelWrapper(model, None, "Escherichia coli",
                                                           'BIOMASS_Ec_iJO1366_core_53p95M')
    return storage._MODELS


@pytest.fixture(scope="function")
def e_coli_core(models):
    """
    Provide the e_coli_core model in a context manager, so that modifications are not persisted beyond the scope of the
    test function. This model is fairly small and should be preferred in test cases where possible.
    """
    with models['test_e_coli_core'].model as model:
        yield model, models['test_e_coli_core'].biomass_reaction


@pytest.fixture(scope="function")
def iJO1366(models):
    """
    Provide the iJO1366 model in a context manager, so that modifications are not persisted beyond the scope of the
    test function.
    """
    with models['test_iJO1366'].model as model:
        yield model, models['test_iJO1366'].biomass_reaction
