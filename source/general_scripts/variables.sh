#!/usr/bin/env bash

BASEDIR="/apps/phenDB"
PHENDB_QUEUE='phenDB'
DB="phenDB"
PHENDB_DEBUG=False

export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"