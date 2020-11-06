#!/usr/bin/env bash

export BASEDIR='/root/rebuild/phenDB'
export PHENDB_QUEUE='phenDB_testing'
export PHENDB_MODEL_DIR='/apps/PICA/models_devel_LL'
export DB='phenDB_devel_LL'
export PHENDB_DEBUG=True

export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
