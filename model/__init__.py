import logging
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
