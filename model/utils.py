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
