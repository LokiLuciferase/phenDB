#!/usr/bin/env bash
set -euo pipefail

cd ${BASEDIR}/source/general_scripts
python3 add_enogs_to_db.py
