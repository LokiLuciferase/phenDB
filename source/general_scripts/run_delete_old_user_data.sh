#!/usr/bin/env bash

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

#export PYTHONPATH="/apps/phenDB_devel_PP/phenDB/source/web_server:$PYTHONPATH"
export PYTHONPATH="/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
# export PYTHONPATH="/apps/phenDB/source/web_server:$PYTHONPATH"

cd /apps/phenDB_devel_LL/source/web_server/businessLogic
# cd /apps/phenDB/source/web_server/businessLogic

python3 delete_old_user_data.py
