FROM biosustain/cameo-solvers:e9e330ca3e37

ENV PYTHONUNBUFFERED=1

ARG CWD=/app
WORKDIR "${CWD}/"

ENV PYTHONPATH=${CWD}/src

RUN pip install --upgrade pip setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade -r /tmp/requirements.txt

COPY . "${CWD}/"

CMD ["gunicorn", "-c", "gunicorn.py", "model.wsgi:app""]
