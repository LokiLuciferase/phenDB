#!/usr/bin/env bash

# cleans database of sample specific data among test runs while leaving enogs and models untouched
#mysql -u root -e "use phenDB; delete from phenotypePredictionApp_bins_in_uploadedfile;"
#mysql -u root -e "use phenDB; delete from phenotypePredictionApp_result_enog;"
#mysql -u root -e "use phenDB; delete from phenotypePredictionApp_result_model;"
#mysql -u root -e "use phenDB; delete from phenotypePredictionApp_bin;"
#mysql -u root -e "use phenDB; delete from phenotypePredictionApp_uploadedfile;"

mysql -u root -e "use phenDB_devel_PP; delete from phenotypePredictionApp_bins_in_uploadedfile;"
mysql -u root -e "use phenDB_devel_PP; delete from phenotypePredictionApp_result_enog;"
mysql -u root -e "use phenDB_devel_PP; delete from phenotypePredictionApp_result_model;"
mysql -u root -e "use phenDB_devel_PP; delete from phenotypePredictionApp_bin;"
mysql -u root -e "use phenDB_devel_PP; delete from phenotypePredictionApp_uploadedfile;"
