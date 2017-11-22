#!/usr/bin/env bash

cd /apps/phenDB/source/web_server/phenotypePredictionApp

sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_UploadedFile;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_model;"
sqlite3 phenDB.sqlite3 "delete from phenotypePredictionApp_result_enog;"
sqlite3 phenDB.sqlite3 "VACUUM;"
