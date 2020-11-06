#!/usr/bin/env bash
set -euo pipefail

cd ${BASEDIR}/source/general_scripts
python3 add_taxonomy_to_db.py
