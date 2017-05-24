FROM biosustain/cameo-solvers:647d6ebdf3dd

COPY requirements.txt requirements.txt
RUN pip install --upgrade --process-dependency-links -r requirements.txt

RUN pip freeze

COPY . ./model
WORKDIR model

ENV PYTHONPATH $PYTHONPATH:/model

ENTRYPOINT ["gunicorn"]
CMD ["-w", "4", "-b", "0.0.0.0:8000", "-t", "150", "-k", "aiohttp.worker.GunicornWebWorker", "model.app:app"]
