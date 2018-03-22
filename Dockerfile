FROM biosustain/cameo-solvers:e9e330ca3e37

# pin pip v9.0.1 to avoid this issue: https://github.com/pypa/pip/issues/5079
# This cobrapy PR should fix the issue: https://github.com/opencobra/cobrapy/pull/680
# When merged, the v9.0.1 pin should be removed.
RUN pip install --upgrade pip==9.0.1 setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade --process-dependency-links -r /tmp/requirements.txt
RUN pip freeze
ADD . /model
WORKDIR /model

ENV PYTHONPATH "${PYTHONPATH}:/model"

CMD ["gunicorn", "-w", "4", "-b", "0.0.0.0:8000", "-t", "150", "-k", "aiohttp.worker.GunicornWebWorker", "model.app:get_app()"]
