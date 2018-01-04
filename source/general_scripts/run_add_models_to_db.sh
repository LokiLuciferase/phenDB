#!/usr/bin/env bash
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="/scratch/swe_ws17/phenDB/source/web_server:$PYTHONPATH"
source activate py3env
cd /scratch/swe_ws17/phenDB/source/general_scripts

python3 add_models_to_db.py