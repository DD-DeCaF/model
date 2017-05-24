FROM biosustain/cameo-solvers:647d6ebdf3dd

RUN pip install --upgrade pip setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade --process-dependency-links -r /tmp/requirements.txt

COPY . /opt/model
WORKDIR /opt/model

ENV PYTHONPATH "${PYTHONPATH}:/opt/model"

ENTRYPOINT ["gunicorn"]
CMD ["-w", "4", "-b", "0.0.0.0:8000", "-t", "150", "-k", "aiohttp.worker.GunicornWebWorker", "model.app:app"]
