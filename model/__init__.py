import logging
import requests
import sys

from raven import Client
from raven.conf import setup_logging
from raven.handlers.logging import SentryHandler

from . import settings

logger = logging.getLogger(__name__)
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
logger.setLevel(logging.INFO)

# Configure Raven to capture warning logs
raven_client = Client(settings.SENTRY_DSN)
handler = SentryHandler(raven_client)
handler.setLevel(logging.WARNING)
setup_logging(handler)

# Retrieve JWT public key from authentication service
try:
    response = requests.get(settings.JWT_PUBLIC_KEY_URL)
    response.raise_for_status()
    settings.JWT_PUBLIC_KEY = response.body
except requests.exceptions.RequestException:
    logger.critical("Could not retrieve JWT public key from authentication service",
                    exc_info=sys.exc_info())
    raise
