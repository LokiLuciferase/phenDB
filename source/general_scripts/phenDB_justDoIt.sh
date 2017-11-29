#!/usr/bin/env bash

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
    cd /apps/phenDB/source/web_server
    nohup chromium-browser http://127.0.0.1:8000/phendb &> /dev/null &
    gunicorn phenotypePrediction.wsgi
}

function showbins() {
    echo "Viewing table 'bin' from database, then exiting..."
    sqlite3 /apps/phenDB/source/web_server/phenotypePredictionApp/phenDB.sqlite3 "select * from phenotypePredictionApp_bin;"
    echo ""
}

function showjobs() {
    echo "Viewing table 'UploadedFile' from database, then exiting..."
    sqlite3 /apps/phenDB/source/web_server/phenotypePredictionApp/phenDB.sqlite3 "select * from phenotypePredictionApp_UploadedFile;"
    echo ""
}

function purge() {
    echo "Purging samples from database..."
    bash /apps/phenDB/source/general_scripts/purge_samples_from_db.sh
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
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done