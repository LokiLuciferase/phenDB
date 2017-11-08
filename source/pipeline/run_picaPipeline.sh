#!/usr/bin/env bash
set -e

CUTOFF="0.75"
INFOLDER="/scratch/swe_ws17/data/test"
JOBNAME=$(basename $INFOLDER)
ABOVE_WORKFOLDER="/scratch/swe_ws17/phenDB_lueftinger/results/"  # for testing, replace this with your testing folder!
WORKFOLDER="${ABOVE_WORKFOLDER}/${JOBNAME}_results"
LOGFOLDER="$WORKFOLDER/logs"
PROGRESS="$LOGFOLDER/progress.log"
LOGLOC="$LOGFOLDER/nextflow.log"
FASTAFILECOUNTFOLDER="$LOGFOLDER/fastafilecount.log"
NODEOFFS=""

mkdir -p $LOGFOLDER
touch $LOGLOC
touch $PROGRESS
# delete also here just to be sure
> $FASTAFILECOUNTFOLDER

module unload java
module load java/1.8u152
module load nextflow

nohup nextflow ./picaPipeline.nf --accuracy_cutoff $CUTOFF --inputfolder $INFOLDER \
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
