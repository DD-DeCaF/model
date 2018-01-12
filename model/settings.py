import os

import requests


ANNOTATIONS_API = os.environ['ANNOTATIONS_API']
ID_MAPPER_API = os.environ['ID_MAPPER_API']
SENTRY_DSN = os.environ.get('SENTRY_DSN', '')

JWT_PUBLIC_KEY = None  # note: will be assigned on init
JWT_PUBLIC_KEY_URL = os.environ['JWT_PUBLIC_KEY_URL']
