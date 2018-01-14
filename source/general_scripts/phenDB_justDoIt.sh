#!/usr/bin/env bash


#BASEDIR="/apps/phenDB"
BASEDIR="/apps/phenDB_devel_LL"
#BASEDIR="/apps/phenDB_devel_PP/phenDB"

#DB="phenDB"
DB="phenDB_devel_LL"
#DB="phenDB_devel_PP"


function usage()
{
    echo "phenDB Utility Script"
    echo ""
    echo -e "\tNo arguments: start chromium browser and run gunicorn server for phenDB."
    echo -e "\t-h --help"
    echo -e "\t--view-bins\tdisplay bins in database"
    echo -e "\t--view-jobs\tdisplay jobs in database"
    echo -e "\t--purge\tremove jobs, bins and associated data\n\t\tfrom database, then start chromium and run server"
    echo ""
}

function runserver() {
    echo "Starting web server and opening a browser window..."
    cd ${BASEDIR}/source/web_server

    if [[ $(pgrep redis-server) = "" ]] || [[ $(ps aux | grep "/usr/bin/[r]q") = "" ]]; then
        echo "Either redis-server or python-rq worker are not running. Exiting."
        exit 1
    fi
    nohup chromium-browser http://127.0.0.1:8000/phendb &> /dev/null &
    python3 manage.py runserver
}

function server_detached() {
    echo "Starting web server as daemon..."
    cd ${BASEDIR}/source/web_server
    if [[ $(pgrep redis-server) = "" ]] || [[ $(ps aux | grep "/usr/bin/[r]q") = "" ]]; then
        echo "Either redis-server or python-rq worker are not running. Exiting."
        exit 1
    fi
    nohup python3 manage.py runserver &>> ${BASEDIR}/logs/django_development_server.log &
}

function showbins() {
    echo "Viewing table 'bin' from database, then exiting..."
    mysql -u root -e "use $DB; select * from phenotypePredictionApp_bin;"
    echo ""
}

function showjobs() {
    echo "Viewing table 'UploadedFile' from database, then exiting..."
    mysql -u root -e "use $DB; select * from phenotypePredictionApp_uploadedfile;"
    echo ""
}

function purge() {
    echo "Purging samples from database..."
    bash ${BASEDIR}/source/general_scripts/purge_samples_from_db.sh
}

function start_queue() {
    echo "Starting queueing tools: redis and python-rq..."
    bash ${BASEDIR}/source/general_scripts/run_redis.sh
}

function monitor_queue() {
    echo "Running the redis monitoring script..."
    python3 ${BASEDIR}/source/general_scripts/monitor_queue.py
}

if [[ $# -eq 0 ]]; then
    runserver
    exit 0
fi

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    case $PARAM in
        -h | --help)
            usage
            exit
            ;;
        --purge)
            purge
            runserver
            ;;
        --view-jobs)
            showjobs
            ;;
        --view-bins)
            showbins
            ;;
        --server-detached)
            server_detached
            ;;
        --start-queue)
            start_queue
            ;;
        --monitor-queue)
            monitor_queue
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done