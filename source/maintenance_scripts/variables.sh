#!/usr/bin/env bash

# these should be fine as is
export BASEDIR="${BASEDIR:-$(pwd)}"
export PHENDB_BASEDIR="${PHENDB_BASEDIR:-$BASEDIR}"
export PYTHONPATH="${PHENDB_BASEDIR}/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export DJANGO_SECRET_KEY="${DJANGO_SECRET_KEY:-testing}"
echo "PHENDB_BASEDIR is now ${PHENDB_BASEDIR}"

# source file of secret variables if exists
SECRETS_FILE="${PHENDB_BASEDIR}/source/maintenance_scripts/secrets.sh"
[ -f "$SECRETS_FILE" ] && source $SECRETS_FILE

# override these if necessary
export PHENDB_DESIGNATION="${PHENDB_DESIGNATION:-phenDB_devel_LL}"
export PHENDB_DATA_DIR="${PHENDB_DATA_DIR:-${PHENDB_BASEDIR}/data}"
export PHENDB_LOG_DIR="${PHENDB_LOG_DIR:-${PHENDB_DATA_DIR}/logs}"
export PHENDB_MODEL_DIR="${PHENDB_MODEL_DIR:-${PHENDB_DATA_DIR}/phenotrex/models}"
export PHENDB_KRONA_DIR="${PHENDB_KRONA_DIR:-${PHENDB_DATA_DIR}/krona/taxonomy}"
export PHENDB_ENOG_ANNOT_FILE="${PHENDB_ENOG_ANNOT_FILE:-${PHENDB_DATA_DIR}/eggnog5/annotations.tsv.gz}"
export PHENDB_ENOG_NAMES_FILE="${PHENDB_ENOG_NAMES_FILE:-${PHENDB_DATA_DIR}/eggnog5/hmmer.db.names.gz}"
export PHENDB_WEB_PORT="${PHENDB_WEB_PORT:-80}"
export PHENDB_DB_PORT="${PHENDB_DB_PORT:-33060}"
export PHENDB_QUEUE="${PHENDB_QUEUE:-${PHENDB_DESIGNATION}}"
export PHENDB_DB_NAME="${PHENDB_DB_NAME:-${PHENDB_DESIGNATION}}"
export PHENDB_DB_USERNAME="${PHENDB_DB_USERNAME:-${PHENDB_DESIGNATION}}"
export PHENDB_DB_PW="${PHENDB_DB_PW:-${PHENDB_DESIGNATION}}"
export PHENDB_ANNOTATION_STRATEGY="${PHENDB_ANNOTATION_STRATEGY:-deepnog}"
export PHENDB_DEBUG="${PHENDB_DEBUG:-True}"
mkdir -p ${PHENDB_LOG_DIR}
