from functools import wraps
import json
import time
import os

from model.logger import logger

def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        ts = time.time()
        result = f(*args, **kw)
        te = time.time()
        logger.info('func:%r args:[%r, %r] took: %2.4f sec' % \
          (f.__name__, args, kw, te-ts))
        return result
    return wrap

def map_reactions_list(map_path):
    """Extract reaction ids from map for FVA optimization

    :param map_path: string
    :return: list of strings
    """
    if not os.path.isfile(map_path):
        return []
    with open(map_path) as f:
        return [i['bigg_id'] for i in json.load(f)[1]['reactions'].values()]
