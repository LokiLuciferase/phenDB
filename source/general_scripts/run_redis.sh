#!/usr/bin/env bash

BASEDIR="/apps/phenDB"
#BASEDIR="/apps/phenDB_devel_LL"
#BASEDIR="/apps/phenDB_devel_PP/phenDB"

mkdir -p ${BASEDIR}/logs
cd ${BASEDIR}/logs

if [[ $(pgrep redis-server) != "" ]] || [[ $(ps aux | grep "/usr/bin/[r]q") != "" ]]; then
    echo "Either redis-server or python-rq worker are already running. Exiting."
    exit 1
fi

# start a redis server and write to log
nohup redis-server &>> redis_server.log &

# start a rq worker and import from the correct location
nohup rq worker --path ${BASEDIR}/source/web_server --path ${BASEDIR}/source/web_server/businessLogic --name phenDB phenDB &>> rq_worker.log &