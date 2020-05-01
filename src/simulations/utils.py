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

import logging
import time
from contextlib import contextmanager
from functools import wraps


logger = logging.getLogger(__name__)


class Singleton(type):
    """Implement the singleton pattern by using this type as metaclass."""

    _instances = {}

    def __call__(cls, *args, **kwargs):
        if cls not in cls._instances:
            cls._instances[cls] = super().__call__(*args, **kwargs)
        return cls._instances[cls]


def timing(f):
    @wraps(f)
    def wrap(*args, **kwargs):
        with log_time(operation=f"func: {f.__name__}, args: [{args}, {kwargs}]"):
            return f(*args, **kwargs)

    return wrap


@contextmanager
def log_time(level=logging.INFO, operation="Task"):
    time_start = time.time()
    yield
    time_end = time.time()
    logger.log(
        level, "{}: completed in {:.4f}s".format(operation, time_end - time_start)
    )
