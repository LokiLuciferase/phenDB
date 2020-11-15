#!/usr/bin/env bash
set -euo pipefail


function usage()
{
    echo "phenDB Utility Script"
    echo ""
    echo -e "\t-h, --help\tDisplay this message and exit"
    echo -e "\t--purge\tremove jobs, bins and associated data\n\t\tfrom database"
    echo -e "\t--start-queue\tActivate redis and python-rq tools if they are not running. Logs at ${PHENDB_LOG_DIR}"
    echo -e "\t--start-debug-server\tActivate Django debug HTTP server for testing. Logs at ${PHENDB_LOG_DIR}"
    echo -e "\t--monitor-queue\tRuns a script to display current status of redis queue."
    echo -e "\t--stop\tStops queue gracefully. Pending jobs are saved."
    echo -e "\t--force-stop\tStops queue immediately. All pending jobs are lost."
    echo ""
}

function purge() {
    echo "Purging samples from database ${DB}..."
    read -p "Continue (Y/N)? " do_continue
    if ! [[ ${do_continue} = "Y" ]]; then
        exit 0
    fi
    python3 ${PHENDB_BASEDIR}/source/maintenance_scripts/purge_user_data.py
}

function start_queue() {
    echo "Starting queueing tools: redis and python-rq..."
    mkdir -p ${PHENDB_LOG_DIR}
    cd ${PHENDB_LOG_DIR}

    if [[ $(pgrep redis-server) != "" ]] || [[ $(ps aux | grep "/usr/bin/[r]q") != "" ]]; then
        echo "Either redis-server or python-rq worker are already running. Exiting."
        exit 1
    fi

    # start a redis server and write to log
    nohup redis-server &>> redis_server.log &

    # wait for redis queue to startup
    sleep 5
    while [[ $(pgrep redis-server) = "" ]]
        do
            sleep 3
        done

    # start a rq worker and import from the correct location; write log to log folder
    nohup rq worker \
    --path ${PHENDB_BASEDIR}/source/web_server \
    --path ${PHENDB_BASEDIR}/source/web_server/businessLogic \
    --name ${PHENDB_QUEUE} ${PHENDB_QUEUE} &>> rq_worker.log &
}


function start_debug_server() {
  echo "Starting debug http server..."
  mkdir -p ${PHENDB_LOG_DIR}
  cd ${PHENDB_LOG_DIR}

  if [[ $(ps aux | grep "[r]unserver") != "" ]]; then
      echo "An instance of the Django development server is already running."
      exit 1
  fi

  nohup python3 ${PHENDB_BASEDIR}/source/web_server/manage.py runserver localhost:8999 &> debug_server.log &

}

function monitor_queue() {
    echo "Running the redis monitoring script..."
    python3 ${PHENDB_BASEDIR}/source/general_scripts/monitor_queue.py
}

# hard shutdown: kill everything
function force_stop() {
    kill $(pgrep -f "manage.py runserver" | tr "\n" " ") || true
    kill $(pgrep redis-server)
}

# soft shutdown: rq finalizes the most current task and saves any queued tasks for later
function stop() {
    echo "Shutting down rq-worker and Redis queue for phenDB..."
    kill $(ps aux | grep usr/bin/[r]q | tr -s " " | cut -f2 -d" ") || true
    while [[ $(ps aux | grep "/usr/bin/[r]q") != "" ]]; do
        sleep 3
    done
    kill $(pgrep -f redis-server | tr "\n" " ") || true
    kill $(pgrep -f "manage.py runserver" | tr "\n" " ") || true
}

# Start redis queue and RQ worker if not yet running.
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
        --start-debug-server)
            start_debug_server
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
