#!/usr/bin/env bash

#export PYTHONPATH=$PYTHONPATH:$HOME/local/lib/python3.6/site-packages

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="/scratch/swe_ws17/phenDB_lueftinger/source/web_server:$PYTHONPATH"

/home/user/lueftinger/miniconda3/envs/py3env/bin/python3 add_models_to_db.py