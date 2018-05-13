#!/usr/bin/env bash

source variables.sh

export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"

cd ${BASEDIR}/source/web_server/businessLogic
python3 precalc_bins.py
