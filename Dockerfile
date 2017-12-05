FROM biosustain/cameo-solvers:e9e330ca3e37

RUN pip install --upgrade pip setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade --process-dependency-links -r /tmp/requirements.txt
RUN pip freeze
ADD . /model
WORKDIR /model

ENV PYTHONPATH "${PYTHONPATH}:/model"

ENTRYPOINT ["gunicorn"]
CMD ["-w", "4", "-b", "0.0.0.0:8000", "-t", "150", "-k", "aiohttp.worker.GunicornWebWorker", "model.app:get_app()"]
