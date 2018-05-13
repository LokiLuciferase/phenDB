#!/usr/bin/env bash

source variables.sh

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

cd ${BASEDIR}/source/web_server/businessLogic
python3 delete_old_user_data.py
