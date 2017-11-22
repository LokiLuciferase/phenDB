#!/usr/bin/env bash
set -e
# TODO: change absolute paths when needed
# check if command line arguments are supplied
if [ $# -eq 3 ]; then
    INFOLDER=$1
    ABOVE_WORKFOLDER=$2
    CUTOFF=$3
else
    INFOLDER="/home/phen_work/singletest"
    ABOVE_WORKFOLDER="/home/phen_work/results/"  # for testing, replace this with your testing folder!
    CUTOFF="0.5"
fi

PIPELINEFILE="/apps/phenDB/source/pipeline/picaPipeline.nf"
export DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
export PYTHONPATH="/apps/phenDB/source/web_server:$PYTHONPATH"

JOBNAME=$(basename $INFOLDER)
WORKFOLDER="${ABOVE_WORKFOLDER}/${JOBNAME}_results"
LOGFOLDER="$WORKFOLDER/logs"
PROGRESS="$LOGFOLDER/progress.log"
LOGLOC="$LOGFOLDER/nextflow.log"
FASTAFILECOUNTFOLDER="$LOGFOLDER/fastafilecount.log"
NODEOFFS=""

mkdir -p $LOGFOLDER
touch $LOGLOC
touch $PROGRESS
# delete also here just to be sure / unnecessary in prod script
> $FASTAFILECOUNTFOLDER

nohup nextflow $PIPELINEFILE --accuracy_cutoff $CUTOFF --inputfolder $INFOLDER \
--workdir $ABOVE_WORKFOLDER --omit_nodes $NODEOFFS -profile standard &> $LOGLOC &

until [ -s $FASTAFILECOUNTFOLDER ]; do
    sleep 1
done

tail -f $LOGLOC