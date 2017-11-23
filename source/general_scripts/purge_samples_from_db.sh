#!/usr/bin/env bash

cd /apps/phenDB/source/web_server/phenotypePredictionApp

# cleans database of sample specific data among test runs while leaving enogs and models untouched
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_UploadedFile;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_bin;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_model;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_enog;"
sqlite3 phenDB.sqlite3 "VACUUM;"
