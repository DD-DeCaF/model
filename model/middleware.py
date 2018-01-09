import logging

from aiohttp import web
import jwt

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


async def auth_middleware(app, handler):
    """rejects the request if a valid token is not provided"""
    async def middleware_handler(request):
        token = request.headers.get('Authorization', '').replace('Bearer ', '')
        try:
            decoded_token = jwt.decode(token, settings.get_jwt_public_key(), algorithms=['RS256'])
        except jwt.exceptions.InvalidTokenError:
            return web.json_response({
                'code': 'unauthorized',
                'message': "Missing or invalid token",
            }, status=401)
        request.jwt_token = decoded_token
        return await handler(request)
    return middleware_handler
