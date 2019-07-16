FROM gcr.io/dd-decaf-cfbf6/modeling-base:master

ENV PYTHONUNBUFFERED=1

ARG CWD=/app
WORKDIR "${CWD}/"

ENV PYTHONPATH=${CWD}/src

# pin pip to 18.0 to avoid issue with cobra -> depinfo -> pipdeptree -> pip._internal.get_installed_distributions
RUN pip install --upgrade pip==18.0 setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade -r /tmp/requirements.txt

COPY . "${CWD}/"
