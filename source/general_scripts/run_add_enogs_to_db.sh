#!/usr/bin/env bash
set -e
source "/apps/phenDB/source/general_scripts/variables.sh"

cd ${BASEDIR}/source/general_scripts
python3 add_enogs_to_db.py
