#!/usr/bin/env bash

set -e
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="/scratch/swe_ws17/phenDB_lueftinger/source/web_server:$PYTHONPATH"

date
echo "Initial size of database: $(du -sh /scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePredictionApp/phenDB.sqlite3)"
/home/user/lueftinger/miniconda3/envs/py3env/bin/python3 /scratch/swe_ws17/phenDB_lueftinger/source/general_scripts/add_enogs_to_db.py

date
echo "Size of db after enog addition: $(du -sh /scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePredictionApp/phenDB.sqlite3)"
/home/user/lueftinger/miniconda3/envs/py3env/bin/python3 /scratch/swe_ws17/phenDB_lueftinger/source/general_scripts/add_models_to_db.py

date
echo "Size of db after model addition: $(du -sh /scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePredictionApp/phenDB.sqlite3)"
echo "Done!"
