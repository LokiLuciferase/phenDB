"""
WSGI config for phenotypePrediction project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/1.11/howto/deployment/wsgi/
"""

import os
import sys

from django.core.wsgi import get_wsgi_application

PHENDB_BASEDIR = os.environ['BASEDIR']
sys.path.append(f"{PHENDB_BASEDIR}/source/web_server")
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "phenotypePrediction.settings")

application = get_wsgi_application()
