#!/usr/bin/env bash

set -e

if [[ "$(hostname)" = "phen.csb.univie.ac.at" ]]; then
    BASEDIR=/apps/phenDB_devel_LL/source
    PY=/usr/bin/python3
else
    BASEDIR=/scratch/swe_ws17/phenDB_lueftinger/source
    PY=/home/user/lueftinger/miniconda3/envs/py3env/bin/python3
fi

function usage()
{
    echo "phenDB Utility Script"
    echo ""
    echo -e "\tNo arguments: start chromium browser and run gunicorn server for phenDB."
    echo -e "\t-d --detached\tall subsequent flags will be processed in detached mode"
    echo -e "\t-h --help"
    echo -e "\t--view-bins\tdisplay bins in database"
    echo -e "\t--view-jobs\tdisplay jobs in database"
    echo -e "\t--purge\tremove jobs, bins and associated data\n\t\tfrom database, then start chromium and run server"
    echo -e "\t--repopulate\tdelete database, re-apply migrations, migrate, repopulate database with models and enogs"
    echo ""
}

function runserver() {
    echo "Starting web server and opening a browser window..."
    cd $BASEDIR/web_server
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

function repopulate() {
    if [[ "$(hostname)" = "phen" ]]; then
        echo "Updating models and enogs is not possible on phen, as there is no eggnog mirror here."
        echo "Please run this script on a vlogin node."
        exit 1
    fi
    echo "Dropping database, and re-populating with models and enogs."
    echo ""
    echo "dropping database..."
    cd $BASEDIR/web_server/phenotypePredictionApp
    rm -rf ./phenDB.sqlite3

    echo "re-applying migrations..."
    cd $BASEDIR/web_server
    $PY manage.py makemigrations
    $PY manage.py migrate

    echo "re-adding enogs..."
    cd $BASEDIR/web_server/phenotypePredictionApp
    export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
    export PYTHONPATH="$BASEDIR/web_server:$PYTHONPATH"

    date
    echo "Initial size of database: $(du -sh ./phenDB.sqlite3)"
    $PY $BASEDIR/general_scripts/add_enogs_to_db.py
    date
    echo "Size of db after enog addition: $(du -sh ./phenDB.sqlite3)"
    echo ""
    echo "re-adding models..."
    $PY $BASEDIR/general_scripts/add_models_to_db.py
    date
    echo "Size of db after model addition: $(du -sh ./phenDB.sqlite3)"
}

if [[ $# -eq 0 ]]; then
    runserver
    exit 0
fi

while [ "$1" != "" ]; do
    PARAM=`echo $1 | awk -F= '{print $1}'`
    case $PARAM in
        -d | --detached)
            shift
            echo "Restarting in detached mode..."
            nohup $BASEDIR/general_scripts/phenDB_justDoIt.sh $* &
            exit 0
            ;;
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
        --repopulate)
            repopulate
            ;;
        *)
            echo "ERROR: unknown parameter \"$PARAM\""
            usage
            exit 1
            ;;
    esac
    shift
done