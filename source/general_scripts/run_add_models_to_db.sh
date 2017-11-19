#!/usr/bin/env bash

#export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python3.6/site-packages

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="/scratch/swe_ws17/phenDB_lueftinger/source/web_server:$PYTHONPATH"
source activate py3env
cd /scratch/swe_ws17/phenDB_lueftinger/source/general_scripts

python3 add_models_to_db.py