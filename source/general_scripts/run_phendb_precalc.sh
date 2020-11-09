#!/usr/bin/env bash
set -euo pipefail

cd ${BASEDIR}/source/web_server/businessLogic
python3 precalc_bins.py
