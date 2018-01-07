#!/usr/bin/env bash

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

#export PYTHONPATH="/apps/phenDB_devel_PP/phenDB/source/web_server:$PYTHONPATH"
#export PYTHONPATH="/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
export PYTHONPATH="/apps/phenDB/source/web_server:$PYTHONPATH"

#source activate py3env
#cd /apps/phenDB_devel_LL/source/general_scripts
#cd /apps/phenDB_devel_PP/phenDB/source/general_scripts
cd /apps/phenDB/source/general_scripts
python3 add_enogs_to_db.py