#!/usr/bin/env python3
import subprocess
import os.path
import threading
from phenotypePredictionApp.models import UploadedFile

class startProcessThread(threading.Thread):
    def __init__(self, keyname):
        threading.Thread.__init__(self)
        self.keyname = keyname
    def run(self):
        print('startProcess called')
        relFilePath = os.path.dirname(UploadedFile.objects.get(key=self.keyname).fileInput.url)
        absPath = os.getcwd()
        infolder = absPath + "/" + relFilePath
        print(infolder)

        #PF: fake script to test webserver
        subprocess.run(["python3", "./fakeScript.py", "--infolder", infolder, "--key", self.keyname])

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