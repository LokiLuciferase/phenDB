#!/usr/bin/env bash

source source/maintenance_scripts/variables.sh
[[ "$(mount | grep /apps/phenDB/data)" == "" && "${PHENDB_DEBUG}" != "True" ]] && echo "Please mount the data directory to /apps/phenDB/data." 1>&2 && return
service mysql restart
bash source/maintenance_scripts/phenDB.sh --load-db
bash source/maintenance_scripts/phenDB.sh --run-fg
bash source/maintenance_scripts/phenDB.sh --dump-db
