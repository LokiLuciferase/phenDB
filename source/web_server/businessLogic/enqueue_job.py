#
# Created by Lukas Lüftinger on 14/01/2018.
#
import sys
import os
import os.path
import shutil
import subprocess
from phenotypePredictionApp.variables import *


def clean_up_on_pipeline_fail(keyname, ppath):

    import django
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    sys.path.append(str(ppath))
    django.setup()
    from phenotypePredictionApp.models import UploadedFile, bins_in_UploadedFile
    currentjob = UploadedFile.objects.get(key=keyname)
    currentjob.errors = None  # set error for website
    # delete bins belonging to failed job
    assoc_rows = bins_in_UploadedFile.objects.filter(UploadedFile=currentjob)
    bins_of_failed = [x.bin for x in assoc_rows]
    for b in bins_of_failed:
        b.delete()
    assoc_rows.delete()
    currentjob.save()


# delete temporary files and uploads after each finished job
def remove_temp_files(infolder=None):

    logfolder = os.path.join(PHENDB_BASEDIR, 'logs')
    shutil.rmtree(os.path.join(logfolder, "work"))
    if infolder:
        shutil.rmtree(infolder)


# delete user submitted data after days
def delete_user_data(days):

    import django
    from django.utils.timezone import make_aware
    from datetime import datetime, timedelta
    from time import time
    import glob

    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    sys.path.append(os.path.join(PHENDB_BASEDIR, "source", "web_server"))
    django.setup()
    from phenotypePredictionApp.models import UploadedFile, bin, result_enog, result_model

    oldest = datetime.today() - timedelta(days=days)

    # delete UploadedFiles older than oldest; association rows deleted automatically
    aged_jobs = UploadedFile.objects.filter(job_date__lte=make_aware(oldest))
    non_precalc_aged = aged_jobs.exclude(key__icontains="PHENDB_PRECALC")
    non_precalc_aged.delete()

    # look for unassociated bins and delete those too
    # spare those bins that have been pre-calculated
    orphan_bins = bin.objects.filter(bins_in_uploadedfile=None)
    orphan_bins.delete()

    # look for unassociated result_enog and result_model
    orphan_hmmer = result_enog.objects.filter(bin=None)
    orphan_verdict = result_model.objects.filter(bin=None)
    orphan_hmmer.delete()
    orphan_verdict.delete()

    # delete results flat files older than days
    oldest_unixtime = float(time()) - (timedelta(days=days).total_seconds())
    result_folders = glob.glob(os.path.join(PHENDB_BASEDIR, "data/results/*"))

    for folder in result_folders:
        if float(os.path.getmtime(folder)) <= float(oldest_unixtime):
            shutil.rmtree(folder)


def phenDB_enqueue(ppath, pipeline_path, infolder, outfolder, pica_cutoff, node_offs):

    # set environmental variables
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    os.environ["PYTHONPATH"] = str(ppath)

    # set pipeline arguments
    arguments = "nextflow {pp} --inputfolder {inf} --outdir {otf} --accuracy_cutoff {pco} \
        --omit_nodes {no} -profile standard".format(pp=pipeline_path,
                                                                 inf=infolder,
                                                                 otf=outfolder,
                                                                 pco=pica_cutoff,
                                                                 no=node_offs)
    # call subprocess and wait for result
    with open(os.path.join(outfolder, "logs/nextflow.log"), "w") as logfile:
        pipeline_call = subprocess.Popen(arguments.split(), stdout=logfile, stderr=logfile)
        pipeline_call.wait()

        # if pipeline encounters error, set errors to None in the DB
        key = os.path.basename(infolder)
        if pipeline_call.returncode != 0:
            remove_temp_files()  # don't delete infolder if failure reason is unknown
            clean_up_on_pipeline_fail(key, ppath)
            raise RuntimeError("Pipeline run has failed. Error status in DB has been updated.")

        remove_temp_files(infolder)

        # if no output file has been generated despite no pipeline error
        # (= if error checking consumes all input files)
        if not os.path.exists(os.path.join(outfolder, "{k}.zip".format(k=key))):
            clean_up_on_pipeline_fail(key, ppath)
            raise RuntimeError("No input files have passed error checking.")

        return pipeline_call.returncode
