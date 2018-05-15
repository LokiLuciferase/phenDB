#!/usr/bin/env bash
set -e
source variables.sh

mkdir -p ${BASEDIR}/logs
cd ${BASEDIR}/logs

# start a redis server and write to log
#FIXME: should already run due to production server
nohup redis-server &>> redis_server.log &

sleep 5
while [[ $(pgrep redis-server) = "" ]]
    do
        sleep 3
    done

# start a rq worker and import from the correct location
nohup rq worker --path ${BASEDIR}/source/web_server --path ${BASEDIR}/source/web_server/businessLogic --name phenDB phenDB &>> rq_worker.log &
