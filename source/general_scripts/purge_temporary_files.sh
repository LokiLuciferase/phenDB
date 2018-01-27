#!/usr/bin/env bash

#BASEDIR="/apps/phenDB"
BASEDIR="/apps/phenDB_devel_LL"

cd ${BASEDIR}/logs

rm -rf work
rm -rf report*
rm -rf traces
rm -rf .nextflow*

cd ../data

rm -rf uploads
