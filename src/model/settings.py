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

import os

import requests


class Settings:
    ENVIRONMENT = os.environ['ENVIRONMENT']
    assert ENVIRONMENT in ('production', 'staging', 'testing', 'development')

    ICE_API = os.environ['ICE_API']
    ICE_USERNAME = os.environ['ICE_USERNAME']
    ICE_PASSWORD = os.environ['ICE_PASSWORD']
    ID_MAPPER_API = os.environ['ID_MAPPER_API']
    MODEL_STORAGE_API = os.environ['MODEL_STORAGE_API']
    SENTRY_DSN = os.environ.get('SENTRY_DSN', '')
    REDIS_ADDR = os.environ['REDIS_ADDR']
    REDIS_PORT = os.environ['REDIS_PORT']
    JWT_PUBLIC_KEY = requests.get(f"{os.environ['IAM_API']}/keys").json()["keys"][0]

    LOGGING = {
        'version': 1,
        'disable_existing_loggers': False,
        'formatters': {
            'simple': {
                'format': "%(asctime)s [%(levelname)s] [%(name)s] %(filename)s:%(funcName)s:%(lineno)d | %(message)s",
            },
        },
        'handlers': {
            'console': {
                'level': 'DEBUG',
                'class': 'logging.StreamHandler',
                'formatter': 'simple',
            },
        },
        'loggers': {
            # All loggers will by default use the root logger below (and
            # hence be very verbose). To silence spammy/uninteresting log
            # output, add the loggers here and increase the loglevel.
        },
        'root': {
            'level': 'DEBUG',
            'handlers': ['console'],
        },
    }
