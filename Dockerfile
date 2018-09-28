FROM gcr.io/dd-decaf-cfbf6/modeling-base:master

ENV PYTHONUNBUFFERED=1

ARG CWD=/app
WORKDIR "${CWD}/"

ENV PYTHONPATH=${CWD}/src

RUN pip install --upgrade pip setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade -r /tmp/requirements.txt

COPY . "${CWD}/"

CMD ["gunicorn", "-c", "gunicorn.py", "model.wsgi:app""]
