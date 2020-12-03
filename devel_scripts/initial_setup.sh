#!/usr/bin/env bash
set -euo pipefail

if [ ! -d "${PHENDB_DATA_DIR}"/phenotrex/models ]; then
  wget https://fileshare.csb.univie.ac.at/phenDB/phenDB_data/latest.tar.gz \
  && tar -xvf latest.tar.gz && rm latest.tar.gz && mv latest/* ${PHENDB_DATA_DIR}; fi

rm -f source/web_server/phenotypePredictionApp/migrations/0*.py
sudo service mysql restart || sudo service mariadb restart
sudo mysql -f < devel_scripts/set_up_dev.sql
python3 source/web_server/manage.py makemigrations phenotypePredictionApp
python3 source/web_server/manage.py migrate
python3 source/web_server/manage.py collectstatic

mkdir -p ${PHENDB_DATA_DIR}/{results,uploads,logs}
mkdir -p ${PHENDB_DATA_DIR}/krona/taxonomy

bash source/maintenance_scripts/purge_temporary_files.sh
python3 source/maintenance_scripts/purge_user_data.py
python3 source/maintenance_scripts/add_taxonomy_to_db.py ${PHENDB_DATA_DIR}/krona/taxonomy
python3 source/maintenance_scripts/add_enogs_to_db.py ${PHENDB_ENOG_ANNOT_FILE}
python3 source/maintenance_scripts/add_models_to_db.py ${PHENDB_MODEL_DIR} --desc_file ${PHENDB_BASEDIR}/source/pipeline/trait_descriptions.txt

# add precalculated refseq data
PRECALC_DIR=${PHENDB_DATA_DIR}/phenotrex/refseq_precalculated_condensed/
python3 devel_scripts/batch_load_refseq_precalc/02_load_to_db.py \
    --preds ${PRECALC_DIR}/predictions.feather \
    --expl ${PRECALC_DIR}/explanations.feather \
    --ccs ${PRECALC_DIR}/compleconta.feather
