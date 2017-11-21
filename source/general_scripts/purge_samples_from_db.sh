#!/usr/bin/env bash

cd /scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePredictionApp

sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_bin;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_model;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_enog;"
sqlite3 phenDB.sqlite3 "VACUUM;"
