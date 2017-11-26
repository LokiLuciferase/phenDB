#!/usr/bin/env python3
import subprocess
import os
import os.path
import threading
from phenotypePredictionApp.models import UploadedFile

class startProcessThread(threading.Thread):
    def __init__(self, keyname):
        threading.Thread.__init__(self)
        self.keyname = keyname

    def run(self):
        print('startProcess called')
        #web_server_folder = "/apps/phenDB/source/web_server"
        web_server_folder = "/apps/PhenDB_PP_devel/phenDB/source/web_server"
        relFilePath = os.path.dirname(UploadedFile.objects.get(key=self.keyname).fileInput.url)
        infolder = os.path.join(str(web_server_folder), str(relFilePath)[1:])
        print("infolder:", infolder)

        # runscript_path = "/apps/phenDB/source/pipeline/run_picaPipeline.sh"
        # pipeline_path = "/apps/phenDB/source/pipeline/picaPipeline.nf"
        # above_workfolder = "/apps/phenDB/source/web_server/results/resultFiles"
        runscript_path = "/apps/PhenDB_PP_devel/phenDB/source/pipeline/run_picaPipeline.sh"
        pipeline_path = "/apps/PhenDB_PP_devel/phenDB/source/pipeline/picaPipeline.nf"
        above_workfolder = "/apps/PhenDB_PP_devel/phenDB/source/web_server/results/resultFiles"

        
        pica_cutoff = "0.5"
        node_offs = ""

        os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
        #os.environ["PYTHONPATH"] = "/apps/phenDB/source/web_server:$PYTHONPATH"
        os.environ["PYTHONPATH"] = "/apps/PhenDB_PP_devel/phenDB/source/web_server:$PYTHONPATH"

        # create workfolder
        outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn=self.keyname))
        os.makedirs(outfolder, exist_ok=True)

        # create log folder
        logfolder = os.path.join(outfolder, "logs")
        os.makedirs(logfolder)

        # call minimal runscript which performs nohup and output rerouting
        subprocess.run([runscript_path,
                        pipeline_path,
                        infolder,
                        outfolder,
                        pica_cutoff,
                        node_offs], check=True)