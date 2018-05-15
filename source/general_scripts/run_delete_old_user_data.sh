#!/usr/bin/env bash
#set -e
#source variables.sh
BASEDIR="/apps/phenDB"

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="${BASEDIR}/source/web_server:$PYTHONPATH"

cd ${BASEDIR}/source/web_server/businessLogic
python3 delete_old_user_data.py
