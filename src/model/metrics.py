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

import prometheus_client


# REQUEST_TIME: Time spent waiting for outgoing API request to internal or external services
# labels:
#   service: The current service (always 'model')
#   environment: The current runtime environment ('production' or 'staging')
#   endpoint: The path to the requested endpoint without query parameters (e.g. '/models/iJO1366')
REQUEST_TIME = prometheus_client.Histogram(
    'decaf_request_handler_duration_seconds',
    "Time spent in request",
    ['service', 'environment', 'endpoint'])


# API_REQUESTS: Time spent waiting for outgoing API request to internal or external services
# labels:
#   service: The current service (always 'model')
#   environment: The current runtime environment ('production' or 'staging')
#   api_name: A short name for the API service (e.g. 'gene-to-reactions')
#   endpoint: The full URL to the external API (e.g.
#             'http://gene-to-reactions/annotations/genes')
API_REQUESTS = prometheus_client.Histogram(
    'decaf_api_request_duration_seconds',
    "Time spent waiting for outgoing API request to internal or external services",
    ['service', 'environment', 'api_name', 'endpoint'])
