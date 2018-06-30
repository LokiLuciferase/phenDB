#!/usr/bin/env bash

#BASEDIR="/apps/phenDB_devel_LL"
BASEDIR="/apps/phenDB"

DB="phenDB"
#DB="phenDB_devel_new"
#DB="phenDB_devel_LL"

# export PYTHONPATH="/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
export PYTHONPATH="/apps/phenDB/source/web_server:$PYTHONPATH"

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"