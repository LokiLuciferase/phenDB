#!/usr/bin/env bash
set -e
source "/apps/phenDB_devel_LL/source/general_scripts/variables.sh"

cd ${BASEDIR}/source/general_scripts
python3 add_models_to_db.py
