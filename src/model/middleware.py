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

from . import settings
from .metrics import REQUEST_TIME


logger = logging.getLogger(__name__)


async def metrics_middleware(app, handler):
    """Log and time all initiated requests"""
    async def middleware_handler(request):
        logger.debug(f"Handling request: {request.url.relative()}")
        with REQUEST_TIME.labels('model', settings.ENVIRONMENT, request.url.path).time():
            return await handler(request)
    return middleware_handler
