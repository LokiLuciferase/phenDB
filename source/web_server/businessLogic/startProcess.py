#!/usr/bin/env python3
import subprocess
#import os.path
#from phenotypePredictionApp.models import UploadedFile


def startProcess(keyname):

    #relFilePath = UploadedFile.objects.get(key = keyname).fileInput.url
    #absPath = PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
    ## start pipeline runscript with path to input folder and output superdirectory
    runscript_path = "/scratch/swe_ws17/phenDB_lueftinger/source/pipeline/run_picaPipeline.sh"
    infolder = "/scratch/swe_ws17/data/test_working"
    above_workfolder = "/scratch/swe_ws17/phenDB_lueftinger/test_runs/results"
    pica_cutoff = "0.5"

    subprocess.run([runscript_path,
                    infolder,
                    above_workfolder,
                    pica_cutoff])


    #TODO: file watcher for progress information + saving in database of certain files
    pass # remove this

if __name__ == "__main__":

    # execute for testing purposes without server
    startProcess("test")