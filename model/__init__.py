import logging
import sys

logger = logging.getLogger('service-model')
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
logger.setLevel(logging.DEBUG)
