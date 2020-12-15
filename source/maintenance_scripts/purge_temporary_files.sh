#!/usr/bin/env bash
set -euo pipefail

rm -rf ${PHENDB_DATA_DIR}/uploads/*
cd ${PHENDB_DATA_DIR}/logs

rm -rf work
rm -rf report*
rm -rf traces
rm -rf .nextflow*

