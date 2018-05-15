#!/usr/bin/env bash
#set -e
BASEDIR="/apps/PhenDB"

mkdir -p ${BASEDIR}/logs
cd ${BASEDIR}/logs

if [[ $(pgrep redis-server) != "" ]] || [[ $(ps aux | grep "/usr/bin/[r]q") != "" ]]; then
    echo "Either redis-server or python-rq worker are already running. Exiting."
    exit 1
fi

# start a redis server and write to log
nohup redis-server &>> redis_server.log &

sleep 5
while [[ $(pgrep redis-server) = "" ]]
    do
        sleep 3
    done

# start a rq worker and import from the correct location
nohup rq worker --path ${BASEDIR}/source/web_server --path ${BASEDIR}/source/web_server/businessLogic --name phenDB phenDB &>> rq_worker.log &