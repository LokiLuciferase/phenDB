#!/usr/bin/env bash

export BASEDIR='/root/rebuild/phenDB'
export PHENDB_QUEUE='phenDB_testing'
export PHENDB_MODEL_DIR='/apps/phenotrex/models/models_live'
export DB='phenDB_devel_LL'
export PHENDB_DB_NAME='phenDB_devel_LL'
export PHENDB_DB_USERNAME='devel_LL'
export PHENDB_DEBUG=True

export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

mkdir -p "$BASEDIR/logs"