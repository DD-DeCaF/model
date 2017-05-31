FROM biosustain/cameo-solvers:e9e330ca3e37

RUN pip install --upgrade pip setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade --process-dependency-links -r /tmp/requirements.txt

COPY . /opt/model
WORKDIR /opt/model

ENV PYTHONPATH "${PYTHONPATH}:/opt/model"

ENTRYPOINT ["gunicorn"]
CMD ["-w", "4", "-b", "0.0.0.0:8000", "-t", "150", "-k", "aiohttp.worker.GunicornWebWorker", "model.app:app"]
