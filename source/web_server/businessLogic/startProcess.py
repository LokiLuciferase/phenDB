#!/usr/bin/env python3
import subprocess
import os.path
from phenotypePredictionApp.models import UploadedFile


def startProcess(keyname):
    print('startProcess called')
    relFilePath = UploadedFile.objects.get(key=keyname).fileInput.url
    absPath = os.getcwd()
    infolder = absPath + "/" + relFilePath
    print(infolder)

    #uncomment when using on the virtual machine
    '''
    runscript_path = "/scratch/swe_ws17/phenDB_lueftinger/source/pipeline/run_picaPipeline.sh"
    above_workfolder = "/scratch/swe_ws17/phenDB_lueftinger/test_runs/results"
    pica_cutoff = "0.5"

    subprocess.run([runscript_path,
                    infolder,
                    above_workfolder,
                    pica_cutoff])

    # TODO: file watcher for progress information + saving in database of certain files
    # TODO: threaded pipeline call?
    '''