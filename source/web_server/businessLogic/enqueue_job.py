#
# Created by Lukas Lüftinger on 14/01/2018.
#
import sys
import os
import os.path
import subprocess

def clean_up_on_pipeline_fail(keyname):

    import django
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    sys.path.append(str(ppath))
    django.setup()
    from phenotypePredictionApp.models import UploadedFile, bins_in_UploadedFile
    currentjob = UploadedFile.objects.get(key=keyname)
    currentjob.errors = None # set error for website
    bins_in_UploadedFile.objects.filter(UploadedFile=currentjob).delete() # delete bins belonging to failed job
    currentjob.save()

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

        # if pipeline encounters error, set errors to None in the DB
        if pipeline_call.returncode != 0:
            clean_up_on_pipeline_fail(os.path.basename(infolder))
            raise RuntimeError("Pipeline run has failed. Error status in DB has been updated.")

        return pipeline_call.returncode
