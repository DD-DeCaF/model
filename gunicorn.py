# Copyright (c) 2018, Novo Nordisk Foundation Center for Biosustainability,
# Technical University of Denmark.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Configure the gunicorn server."""

import os

from prometheus_client import multiprocess


_config = os.environ["ENVIRONMENT"]

bind = "0.0.0.0:8000"
# This service, as an exception, uses synchronous workers to avoid having to
# share model state between workers. This allows us to use cobrapy's context
# manager while making modifications to the shared model instances, resetting
# them on completion.
worker_class = "sync"
# The timeout is increased from the default of 20s, to allow time to deserialize
# proprietary models that are not preloaded on service startup.
timeout = 120


def child_exit(server, worker):
    multiprocess.mark_process_dead(worker.pid)


if _config in ['production', 'staging']:
    # Considerations when choosing the number of synchronous workers:
    # - How many simultaneous requests do we want to be able to handle
    # - Separate workers will cache non-public models separately, meaning the
    #   more workers the less likely to get a cache hit.
    workers = 6
    preload_app = True
elif _config in ['testing', 'development']:
    workers = 1
    reload = True
    accesslog = "-"
    access_log_format = '''%(t)s "%(r)s" %(s)s %(b)s %(L)s "%(f)s"'''
