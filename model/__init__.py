import logging
import sys

logger = logging.getLogger()
logger.addHandler(logging.StreamHandler(stream=sys.stdout))
logger.setLevel(logging.DEBUG)
