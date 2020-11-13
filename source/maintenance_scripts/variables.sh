#!/usr/bin/env bash

# override these if necessary
export PHENDB_MODEL_DIR="${PHENDB_MODEL_DIR:-/apps/phenotrex/models/models_live}"
export PHENDB_QUEUE="${PHENDB_QUEUE:-phenDB_devel_LL}"
export PHENDB_DB_NAME="${PHENDB_DB_NAME:-phenDB_devel_LL}"
export PHENDB_DB_USERNAME="${PHENDB_DB_USERNAME:-phenDB_devel_LL}"
export PHENDB_DB_PW="${PHENDB_DB_PW:-phenDB_devel_LL}"
export PHENDB_ANNOTATION_STRATEGY="${PHENDB_ANNOTATION_STRATEGY:-hmmer}"
export PHENDB_DEBUG="${PHENDB_DEBUG:-True}"

# these should be fine as is
SCRIPTPATH="$( cd "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
export BASEDIR="$(cd ${SCRIPTPATH}; cd ../..; pwd)"
export PHENDB_BASEDIR=${BASEDIR}
export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

mkdir -p "$BASEDIR/logs"
