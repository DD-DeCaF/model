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
import time

from flask import g, request

from . import settings
from .metrics import REQUEST_TIME


logger = logging.getLogger(__name__)


def init_app(app):
    @app.before_request
    def before_request():
        logger.debug(f"Handling request: {request.path}")
        g.request_start = time.time()

    @app.after_request
    def after_request(response):
        request_duration = time.time() - g.request_start
        REQUEST_TIME.labels('model',
                            app.config['ENVIRONMENT'],
                            request.path).observe(request_duration)
        return response
