#!/usr/bin/env bash
set -e
source "/apps/phenDB/source/general_scripts/variables.sh"

cd ${BASEDIR}/source/web_server/businessLogic
python3 delete_old_user_data.py
