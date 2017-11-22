from functools import wraps
import inspect
import logging
import time


def timing(f):
    @wraps(f)
    def wrap(*args, **kw):
        time_start = time.time()
        result = f(*args, **kw)
        time_end = time.time()
        # Get the appropriate logger for the file

        try:
            logger_name = inspect.getmodule(f).__name__
        except TypeError:
            logger_name = 'builtin'
        function_logger = logging.getLogger(logger_name)
        function_logger.info('func:%r args:[%r, %r] took: %2.4f sec',
                             f.__name__, args, kw, time_end - time_start)
        return result
    return wrap
