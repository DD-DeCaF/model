FROM biosustain/cameo-solvers:e9e330ca3e37

ENV PYTHONUNBUFFERED=1

ARG CWD=/app
WORKDIR "${CWD}/"

ENV PYTHONPATH=${CWD}/src

# pin pip v9.0.1 to avoid this issue: https://github.com/pypa/pip/issues/5079
# This cobrapy PR should fix the issue: https://github.com/opencobra/cobrapy/pull/680
# When merged, the v9.0.1 pin should be removed.
RUN pip install --upgrade pip==9.0.1 setuptools wheel
COPY requirements.txt /tmp/requirements.txt
RUN pip install --upgrade --process-dependency-links -r /tmp/requirements.txt
RUN pip freeze

# FIXME: This ad-hoc patch is necessary to avoid optlang overriding our root logger level.
# Hopefully, this will be included in an optlang release.
# See https://github.com/biosustain/optlang/pull/162
RUN sed -i "s/getLogger()/getLogger('cplex.problem')/" /usr/local/lib/python3.6/site-packages/optlang/cplex_interface.py

COPY . "${CWD}/"

CMD ["gunicorn", "-c", "gunicorn.py", "model.wsgi:app""]
