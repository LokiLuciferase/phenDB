#!/usr/bin/env bash

source "/apps/phenDB_devel_LL/source/general_scripts/variables.sh"

function usage()
{
    echo "phenDB Utility Script"
    echo ""
    echo -e "\tNo arguments: start chromium browser and run django development server for phenDB."
    echo -e "\t-h, --help\tDisplay this message and exit"
    echo -e "\t--view-bins\tdisplay bins in database"
    echo -e "\t--view-jobs\tdisplay jobs in database"
    echo -e "\t--purge\tremove jobs, bins and associated data\n\t\tfrom database"
    echo -e "\t--server-detached\tStart django development server in detached mode - log file at $BASEDIR/logs"
    echo -e "\t--start-queue\tActivate redis and python-rq tools if they are not running. Logs at $BASEDIR/logs"
    echo -e "\t--monitor-queue\tRuns a script to display current status of redis queue."
    echo -e "\t--start\tChecks for running services, then runs start-queue and server-detached."
    echo -e "\t--stop\tStops queue and development web server gracefully. Pending jobs are saved."
    echo -e "\t--force-stop\tStops server and queue immediately. All pending jobs are lost."
    echo ""
}

function purge() {
   echo "Purging samples from database..."
   mysql -u root -e "use ${DB}; delete from phenotypePredictionApp_bininjob;"
   mysql -u root -e "use ${DB}; delete from phenotypePredictionApp_hmmerresult;"
   mysql -u root -e "use ${DB}; delete from phenotypePredictionApp_picaresult;"
   mysql -u root -e "use ${DB}; delete from phenotypePredictionApp_bin;"
   mysql -u root -e "use ${DB}; delete from phenotypePredictionApp_job;"
}

function start_queue() {
   echo "Starting queueing tools: redis and python-rq..."
   bash ${BASEDIR}/source/general_scripts/run_redis.sh
}

function monitor_queue() {
   echo "Running the redis monitoring script..."
   python3 ${BASEDIR}/source/general_scripts/monitor_queue.py
}

function force_stop() {
   # hard shutdown: kill everything
   kill $(ps aux | grep "/usr/bin/python3 manage.py runserve[r]" | tr -s " " | cut -f2 -d" ")
   kill $(pgrep redis-server)
}

function stop() {
   # soft shutdown: rq finalizes the most current task and saves any queued tasks for later
   echo "Shutting down Redis queue and Django Development Server for phenDB..."
   kill $(ps aux | grep "/usr/bin/python3 manage.py runserve[r]" | tr -s " " | cut -f2 -d" ")
   kill $(ps aux | grep usr/bin/[r]q | tr -s " " | cut -f2 -d" ")
   kill $(pgrep redis-server)
}

function start() {
   echo "Checking if service is running..."
   if ! [[ $(pgrep redis-server) = "" ]] || ! [[ $(ps aux | grep "/usr/bin/[r]q") = "" ]]; then
       echo "Either redis-server or python-rq worker are running."
       stop
   fi
   echo "Starting Redis queue..."
   start_queue
   exit 0
}

if [[ $# -eq 0 ]]; then
    usage
    exit 0
fi

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    case ${PARAM} in
        -h | --help)
            usage
            exit
            ;;
        --purge)
            purge
            ;;
        --start-queue)
            start_queue
            ;;
        --monitor-queue)
            monitor_queue
            ;;
        --force-stop)
            force_stop
            ;;
        --start)
            start
            ;;
        --stop)
            stop
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done