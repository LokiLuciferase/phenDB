#
# Created by Lukas Lüftinger on 14/01/2018.
#
import os
import subprocess


def phenDB_enqueue(runscript_path, pipeline_path, infolder, outfolder, pica_cutoff, node_offs):

    # set environmental variables
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    # os.environ["PYTHONPATH"] = "/apps/phenDB/source/web_server:$PYTHONPATH"
    os.environ["PYTHONPATH"] = "/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
    # os.environ["PYTHONPATH"] = "/apps/phenDB_devel_PP/phenDB/source/web_server:$PYTHONPATH"

    # call minimal runscript which performs nohup and output rerouting
    subprocess.run([runscript_path,
                    pipeline_path,
                    infolder,
                    outfolder,
                    pica_cutoff,
                    node_offs], check=True)