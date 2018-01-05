#!/usr/bin/env bash
set -e
LOGLOC="${3}/logs/nextflow.log"
nohup nextflow $1 --inputfolder $2 --outdir $3 --accuracy_cutoff $4 --omit_nodes $5 -profile standard -with-report -with-dag flowchart.png with-timeline &> $LOGLOC &