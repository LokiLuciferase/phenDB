#!/usr/bin/env bash
set -e
source "/apps/phenDB/source/general_scripts/variables.sh"

cd ${BASEDIR}/logs

rm -rf work
rm -rf report*
rm -rf traces
rm -rf .nextflow*

cd ../data

rm -rf uploads
