import subprocess
import os.path
from phenotypePredictionApp.models import UploadedFile


def startProcess(keyname):

    #relFilePath = UploadedFile.objects.get(key = keyname).fileInput.url
    #absPath =  PROJECT_ROOT = os.path.abspath(os.path.dirname(__file__))
    #subprocess.run(["/folder/folder/run_picaPipeline.sh", "--inputfolder", "")


    #TODO: file watcher for progress information + saving in database of certain files
    pass # remove this