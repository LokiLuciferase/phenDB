#
# Created by Lukas Lüftinger on 14/01/2018.
#
import os
import os.path
import subprocess


def phenDB_enqueue(runscript_path, pipeline_path, infolder, outfolder, pica_cutoff, node_offs):

    # set environmental variables
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

    # os.environ["PYTHONPATH"] = "/apps/phenDB/source/web_server:$PYTHONPATH"
    os.environ["PYTHONPATH"] = "/apps/phenDB_devel_LL/source/web_server:$PYTHONPATH"
    # os.environ["PYTHONPATH"] = "/apps/phenDB_devel_PP/phenDB/source/web_server:$PYTHONPATH"

    arguments = "nextflow {pp} --inputfolder {inf} --outdir {otf} --accuracy_cutoff {pco} \
        --omit_nodes {no} -profile standard -with-report".format(pp=pipeline_path,
                                                                 inf=infolder,
                                                                 otf=outfolder,
                                                                 pco=pica_cutoff,
                                                                 no=node_offs)

    with open(os.path.join(outfolder, "logs/nextflow.log"), "w") as logfile:
        pipeline_call = subprocess.Popen(arguments.split(), stdout=logfile, stderr=logfile)
        return pipeline_call.returncode