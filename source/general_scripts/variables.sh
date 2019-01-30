#!/usr/bin/env bash

BASEDIR="/apps/phenDB"
DB="phenDB"

export PYTHONPATH="/apps/phenDB/source/web_server:$PYTHONPATH"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"