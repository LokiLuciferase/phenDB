#
# Created by Lukas Lüftinger on 14/01/2018.
#
import os
import os.path
import subprocess


def phenDB_enqueue(ppath, pipeline_path, infolder, outfolder, pica_cutoff, node_offs):

    # set environmental variables
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    os.environ["PYTHONPATH"] = str(ppath)

    # set pipeline arguments
    arguments = "nextflow {pp} --inputfolder {inf} --outdir {otf} --accuracy_cutoff {pco} \
        --omit_nodes {no} -profile standard -with-report".format(pp=pipeline_path,
                                                                 inf=infolder,
                                                                 otf=outfolder,
                                                                 pco=pica_cutoff,
                                                                 no=node_offs)

    # call subprocess and wait for result
    with open(os.path.join(outfolder, "logs/nextflow.log"), "w") as logfile:
        pipeline_call = subprocess.Popen(arguments.split(), stdout=logfile, stderr=logfile)
        pipeline_call.wait()
        if pipeline_call.returncode != 0:
            raise RuntimeError("Pipeline has exited with error.")
        return pipeline_call.returncode