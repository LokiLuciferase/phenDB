#!/usr/bin/env bash

BASEDIR="/apps/phenDB_devel_LL"
DB="phenDB_devel_LL"

export PYTHONPATH="/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"