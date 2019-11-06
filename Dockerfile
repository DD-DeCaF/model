FROM gcr.io/dd-decaf-cfbf6/modeling-base:master

ARG CWD=/app
WORKDIR "${CWD}/"

ENV PYTHONPATH=${CWD}/src

COPY requirements ./requirements

RUN pip-sync requirements/requirements.txt

COPY . "${CWD}/"
