#!/usr/bin/env bash

export BASEDIR="/apps/phenDB"
export PHENDB_QUEUE='phenDB'
export DB="phenDB"
export PHENDB_DEBUG=False

export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"