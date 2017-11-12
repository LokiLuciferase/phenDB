#!/usr/bin/env bash
set -e
# TODO: change absolute paths when needed
# check if command line arguments are supplied
if [ $# -eq 3 ]; then
    INFOLDER=$1
    ABOVE_WORKFOLDER=$2
    CUTOFF=$3
else
    INFOLDER="/scratch/swe_ws17/data/test"
    ABOVE_WORKFOLDER="/scratch/swe_ws17/phenDB_lueftinger/results/"  # for testing, replace this with your testing folder!
    CUTOFF="0.75"
fi

PIPELINEFILE="/scratch/swe_ws17/phenDB_lueftinger/source/pipeline/picaPipeline.nf"
DJANGO_SETTINGS_MODULE="phenotypePrediction.settings"
PYTHONPATH="/scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePrediction:/scratch/swe_ws17/phenDB_lueftinger/source/web_server/phenotypePredictionApp:$PYTHONPATH"
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

module unload java
module load java/1.8u152
module load nextflow

nohup nextflow $PIPELINEFILE --accuracy_cutoff $CUTOFF --inputfolder $INFOLDER \
--workdir $ABOVE_WORKFOLDER --omit_nodes $NODEOFFS -profile cluster  &> $LOGLOC &
sleep 3

until [ -s $FASTAFILECOUNTFOLDER ]; do
    sleep 1
done
totalnum=$(cat $FASTAFILECOUNTFOLDER)

res="0.0"
while [ "$res" != "100.00" ]; do
    sleep 1
    donenum=$(cat $PROGRESS | wc -l)
    res=$(bc <<< "scale=2; ($donenum/($totalnum))*100")
    echo -ne "Job completion: $res%"\\r
done
