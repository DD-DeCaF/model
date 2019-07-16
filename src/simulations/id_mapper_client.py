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

import json
import logging
import os

import requests

from simulations.app import app
from simulations.metrics import API_REQUESTS
from simulations.utils import log_time


logger = logging.getLogger(__name__)


def query_identifiers(object_ids, db_from, db_to):
    """
    Call the id mapper service.

    :param object_ids: list of identifiers to query
    :param db_from: the source of the identifier, e.g. 'kegg'
    :param db_to: the destination type of the identifier, e.g. 'bigg'
    """
    if len(object_ids) == 0:
        return {}
    query = json.dumps({'ids': object_ids, 'dbFrom': db_from, 'dbTo': db_to, 'type': 'Metabolite'})
    logger.info('query id mapper at %s with %s', app.config['ID_MAPPER_API'], str(query))
    with log_time(operation=f"ID map request for ids: {object_ids}"):
        with API_REQUESTS.labels('model', os.environ['ENVIRONMENT'], 'id-mapper', app.config['ID_MAPPER_API']).time():
            return requests.post(f"{app.config['ID_MAPPER_API']}/query", data=query).json()['ids']
