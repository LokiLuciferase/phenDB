#!/usr/bin/env python3
import subprocess
import os.path
from phenotypePredictionApp.models import UploadedFile


def startProcess(keyname):

    # relFilePath = UploadedFile.objects.get(key = keyname).fileInput.url
    # absPath = PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
    ## start pipeline runscript with path to input folder and output superdirectory
    print('startProcess called')
    relFilePath = UploadedFile.objects.get(key=keyname).fileInput.url
    absPath = os.getcwd()
    infolder = absPath + "/" + relFilePath
    print(infolder)

    runscript_path = "/apps/phenDB/source/pipeline/run_picaPipeline.sh"
    above_workfolder = "/home/phen_work/results"
    pica_cutoff = "0.5"

    subprocess.run([runscript_path,
                    infolder,
                    above_workfolder,
                    pica_cutoff])

    # TODO: file watcher for progress information + saving in database of certain files
    # TODO: threaded pipeline call?
