import logging

from aiohttp import web

from . import settings
from . import raven_client


async def raven_middleware(app, handler):
    """aiohttp middleware captures any uncaught exceptions to Sentry before re-raising"""
    async def middleware_handler(request):
        try:
            return await handler(request)
        except Exception:
            raven_client.captureException()
            raise
    return middleware_handler
