#!/usr/bin/env bash
set -euo pipefail

rm -f source/web_server/phenotypePredictionApp/migrations/0*.py
sudo service mysql start || sudo service mariadb start
sudo mysql < devel_scripts/set_up_dev.sql
python3 source/web_server/manage.py makemigrations phenotypePredictionApp
python3 source/web_server/manage.py migrate

mkdir -p ${PHENDB_DATA_DIR}/{results,uploads,logs}
mkdir -p ${PHENDB_DATA_DIR}/krona/taxonomy

python3 source/maintenance_scripts/purge_user_data.py
python3 source/maintenance_scripts/add_taxonomy_to_db.py ${PHENDB_DATA_DIR}/krona/taxonomy
python3 source/maintenance_scripts/add_enogs_to_db.py ${PHENDB_DATA_DIR}/annotations.tsv.gz
python3 source/maintenance_scripts/add_models_to_db.py ${PHENDB_DATA_DIR}/models
