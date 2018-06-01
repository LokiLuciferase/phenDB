#!/usr/bin/env bash
set -e
source "/apps/phenDB_devel_LL/source/general_scripts/variables.sh"

cd ${BASEDIR}/source/web_server/businessLogic
python3 precalc_bins.py
