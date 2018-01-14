#!/usr/bin/env bash

cd /apps/redis/logs

# start a redis server and write to log
nohup redis-server &>> redis_server.log &

# start a rq worker and import from the correct location
nohup rq worker --path /apps/phenDB/source/web_server --name phenDB phenDB &>> rq_worker.log &