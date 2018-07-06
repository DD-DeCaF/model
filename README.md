# model

![master Branch](https://img.shields.io/badge/branch-master-blue.svg)
[![master Build Status](https://travis-ci.org/DD-DeCaF/model.svg?branch=master)](https://travis-ci.org/DD-DeCaF/model)
[![master Codecov](https://codecov.io/gh/DD-DeCaF/model/branch/master/graph/badge.svg)](https://codecov.io/gh/DD-DeCaF/model/branch/master)
[![master Requirements Status](https://requires.io/github/DD-DeCaF/model/requirements.svg?branch=master)](https://requires.io/github/DD-DeCaF/model/requirements/?branch=master)

![devel Branch](https://img.shields.io/badge/branch-devel-blue.svg)
[![devel Build Status](https://travis-ci.org/DD-DeCaF/model.svg?branch=devel)](https://travis-ci.org/DD-DeCaF/model)
[![devel Codecov](https://codecov.io/gh/DD-DeCaF/model/branch/devel/graph/badge.svg)](https://codecov.io/gh/DD-DeCaF/model/branch/devel)
[![devel Requirements Status](https://requires.io/github/DD-DeCaF/model/requirements.svg?branch=devel)](https://requires.io/github/DD-DeCaF/model/requirements/?branch=devel)

## Development

Type `make` to see frequently used workflow commands.

### Data

To update models, run `make update_models`.

To update salts, run `make update_salts`.

### Testing

To run all tests and QA checks, run `make qa`.

### Environment

Specify environment variables in a `.env` file. See `docker-compose.yml` for the possible variables and their default values.

* `ENVIRONMENT` Set to either `development`, `testing`, `staging` or `production`
* `SENTRY_DSN` DSN for reporting exceptions to [Sentry](https://docs.sentry.io/clients/python/integrations/flask/).
* `REDIS_ADDR` Local Redis service hostname
* `REDIS_PORT` Local Redis service port
* `ICE_API` ICE API endpoint
* `ICE_USERNAME` ICE username
* `ICE_PASSWORD` ICE password
* `ID_MAPPER_API` URL to the ID mapper service
