#!/usr/bin/env bash
set -e
source variables.sh

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

cd ${BASEDIR}/source/general_scripts
python3 add_models_to_db.py
