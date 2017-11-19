import subprocess
import os.path
from phenotypePredictionApp.models import UploadedFile


def startProcess(keyname):
    print('startProcess called')
    relFilePath = UploadedFile.objects.get(key = keyname).fileInput.url
    absPath = os.getcwd()
    path = absPath + "/" + relFilePath
    print(path)
    #subprocess.run(["/folder/folder/run_picaPipeline.sh", "--inputfolder", "")


    #TODO: file watcher for progress information + saving in database of certain files