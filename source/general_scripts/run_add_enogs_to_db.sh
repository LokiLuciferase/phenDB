#!/usr/bin/env bash

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
#export PYTHONPATH="/scratch/swe_ws17/phenDB_lueftinger/source/web_server:$PYTHONPATH"
export PYTHONPATH="/scratch/swe_ws17/phenDB_peneder/source/web_server:$PYTHONPATH"
#source activate py3env
#cd /scratch/swe_ws17/phenDB_lueftinger/source/general_scripts
cd /apps/phenDB_devel_PP/phenDB/source/general_scripts



python3 add_enogs_to_db.py