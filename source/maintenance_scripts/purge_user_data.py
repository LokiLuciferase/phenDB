import os
import sys
import django
from django.utils.timezone import make_aware
from datetime import datetime, timedelta
from time import time
import glob
import shutil

import click

PHENDB_BASEDIR = os.environ['PHENDB_BASEDIR']
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
sys.path.append(os.path.join(PHENDB_BASEDIR, "source", "web_server"))

django.setup()
from phenotypePredictionApp.models import Job, Bin, BinInJob, HmmerResult, PicaResult


@click.command()
@click.option('--days', type=int, default=-1)
def purge_user_data(days):
    oldest = datetime.today() - timedelta(days=days)

    # delete Jobs older than oldest
    # spare those bins and results that have been pre-calculated from refseq genomes
    aged_jobs = Job.objects.filter(job_date__lte=make_aware(oldest))
    non_precalc_aged = aged_jobs.exclude(key__icontains="PHENDB_PRECALC")
    non_precalc_aged.delete()

    # look for unassociated bins and bij and delete those too
    orphan_bij = BinInJob.objects.filter(job=None)
    orphan_bij.delete()

    orphan_bins = Bin.objects.filter(bininjob=None).exclude(bininjob__job__key__icontains="PHENDB_PRECALC")
    orphan_bins.delete()

    # look for unassociated result_enog and result_model
    orphan_hmmer = HmmerResult.objects.filter(bin=None)
    orphan_verdict = PicaResult.objects.filter(bin=None)
    orphan_hmmer.delete()
    orphan_verdict.delete()

    # delete results flat files older than days
    oldest_unixtime = float(time()) - (timedelta(days=days).total_seconds())
    result_basedir = os.path.join(PHENDB_BASEDIR, "data/results/*")
    result_folders = glob.glob(result_basedir)

    for folder in result_folders:
        if float(os.path.getmtime(folder)) <= float(oldest_unixtime):
            shutil.rmtree(folder)

    # delete any file for which no job is registered in the database
    uploadfolder = os.path.join(PHENDB_BASEDIR, "data/uploads")
    retained_upload_folders = os.listdir(uploadfolder)
    for key in retained_upload_folders:
        try:
            retained_upload_job = Job.objects.get(key=key)
        except:
            shutil.rmtree(os.path.join(uploadfolder, key))


if __name__ == '__main__':
    purge_user_data()
