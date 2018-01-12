import asyncio
import aiohttp_cors
from aiohttp import web
import logging

import model.handlers as handlers
from .middleware import raven_middleware, auth_middleware


LOGGER = logging.getLogger(__name__)


def get_app():
    app = web.Application(middlewares=[raven_middleware, auth_middleware])
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
    server = await loop.create_server(app.make_handler(), '0.0.0.0', 8000)
    LOGGER.info('Web server is up')
    return app


if __name__ == '__main__':
    loop = asyncio.get_event_loop()
    loop.run_until_complete(start(loop))
    try:
        loop.run_forever()
    except KeyboardInterrupt:
        pass
