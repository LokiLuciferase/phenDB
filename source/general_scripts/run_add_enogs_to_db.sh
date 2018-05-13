#!/usr/bin/env bash

source variables.sh

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

cd ${BASEDIR}/source/general_scripts
python3 add_enogs_to_db.py
