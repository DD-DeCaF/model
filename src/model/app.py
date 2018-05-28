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

import asyncio
import logging

import aiohttp_cors
from aiohttp import web

import model.handlers as handlers

from .middleware import raven_middleware


logger = logging.getLogger(__name__)


def get_app():
    app = web.Application(middlewares=[raven_middleware])
    app.router.add_route('GET', '/wsmodels/{model_id}', handlers.model_ws_full)
    app.router.add_route('GET', '/v1/wsmodels/{model_id}', handlers.model_ws_json_diff)
    app.router.add_route('GET', '/maps', handlers.maps)
    app.router.add_route('GET', '/map', handlers.map)
    app.router.add_route('GET', '/model-options/{species}', handlers.model_options)
    app.router.add_route('POST', '/models/{model_id}', handlers.model)
    app.router.add_route('GET', '/v1/models/{model_id}', handlers.model_get)
    app.router.add_route('POST', '/v1/models/{model_id}', handlers.model_diff)
    app.router.add_route('GET', '/v1/model-info/{model_id}', handlers.model_info)

    # Configure default CORS settings.
    cors = aiohttp_cors.setup(app, defaults={
        "*": aiohttp_cors.ResourceOptions(
            allow_credentials=True,
            expose_headers="*",
            allow_headers="*",
        )
    })

    # Configure CORS on all routes.
    for route in list(app.router.routes()):
        cors.add(route)
    return app


async def start(loop):
    app = get_app()
    await loop.create_server(app.make_handler(), '0.0.0.0', 8000)
    logger.info('Web server is up')
    return app


if __name__ == '__main__':
    loop = asyncio.get_event_loop()
    loop.run_until_complete(start(loop))
    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass
