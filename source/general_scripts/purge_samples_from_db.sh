#!/usr/bin/env bash

if [[ $# -eq 0 ]]; then
    echo "Please provide database name as argument."
    exit 1
fi

# cleans database of sample specific data among test runs while leaving enogs and models untouched
mysql -u root -e "use $1; delete from phenotypePredictionApp_bins_in_uploadedfile;"
mysql -u root -e "use $1; delete from phenotypePredictionApp_result_enog;"
mysql -u root -e "use $1; delete from phenotypePredictionApp_result_model;"
mysql -u root -e "use $1; delete from phenotypePredictionApp_bin;"
mysql -u root -e "use $1; delete from phenotypePredictionApp_uploadedfile;"
