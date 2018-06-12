#
# Created by Lukas Lüftinger on 14/01/2018.
#
import sys
import os
import os.path
import shutil
import subprocess
from phenotypePredictionApp.variables import PHENDB_BASEDIR, PHENDB_QUEUE, PHENDB_DEBUG


TAXONOMY_NAMES_FILE = "/var/www/krona_taxonomy/taxonomy.tab"  # using kronatools taxonomy names, easily updatable


def update_taxonomy(ppath):
    import django
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    sys.path.append(ppath)
    django.setup()
    from phenotypePredictionApp.models import Taxon

    print("Updating KronaTools Taxonomy DB using updateTaxonomy.sh...")
    update_call = subprocess.run("/apps/KronaTools/updateTaxonomy.sh", check=True)
    if not update_call.returncode == 0:
        raise RuntimeError("Updating of taxonomy.tab exited with error code. Aborting.")
    print("Done.")

    if not os.path.exists(TAXONOMY_NAMES_FILE):
        raise RuntimeError("taxonomy.tab was not found. Cannot update Taxonomy.")

    print("Reading rows from updated taxonomy file...")
    taxonomy_entries = []
    counter = 0
    with open(TAXONOMY_NAMES_FILE, "r") as taxfile:
        for line in taxfile:
            sys.stdout.write("Reading line {i}               \r".format(i=counter))
            sys.stdout.flush()
            counter += 1
            ix, parent, tax_id, rank, namestring = line.strip().split("\t")
            new_taxon = Taxon(tax_id=tax_id, taxon_name=namestring)
            taxonomy_entries.append(new_taxon)
    print("Done.")

    if len(taxonomy_entries) < 1700000:
        raise RuntimeError("Something went wrong during database read-in. Aborting.")


    # drop all rows from old taxonomy DB
    print("Rebuilding taxonomy DB...")
    Taxon.objects.all().delete()
    while len(taxonomy_entries) > 0:
        sys.stdout.write("                {n}                entries left to add.\r".format(n=len(taxonomy_entries)))
        sys.stdout.flush()
        subset = taxonomy_entries[-10000:]
        taxonomy_entries = taxonomy_entries[:-10000]
        Taxon.objects.bulk_create(subset)
    print("Finished updating Taxonomy Table.")


def clean_up_on_pipeline_fail(keyname, ppath, failtype):

    import django
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    sys.path.append(str(ppath))
    django.setup()
    from phenotypePredictionApp.models import Job, Bin, BinInJob
    currentjob = Job.objects.get(key=keyname)
    currentjob.errors = None  # set error for website
    currentjob.error_type = failtype
    # delete bins belonging to failed job
    assoc_rows = BinInJob.objects.filter(job=currentjob)
    bins_of_failed = [x.bin for x in assoc_rows]
    for b in bins_of_failed:
        b.delete()
    assoc_rows.delete()
    currentjob.save()


# delete temporary files and uploads after each finished job
def remove_temp_files(infolder=None):

    if PHENDB_DEBUG:
        return
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
    from phenotypePredictionApp.models import Job, Bin, HmmerResult, PicaResult

    oldest = datetime.today() - timedelta(days=days)

    # delete Jobs older than oldest; association rows deleted automatically
    # spare those Jobs that have been pre-calculated
    aged_jobs = Job.objects.filter(job_date__lte=make_aware(oldest))
    non_precalc_aged = aged_jobs.exclude(key__icontains="PHENDB_PRECALC")
    non_precalc_aged.delete()

    # look for unassociated bins and delete those too
    orphan_bins = Bin.objects.filter(BinInJob=None)
    orphan_bins.delete()

    # look for unassociated result_enog and result_model
    orphan_hmmer = HmmerResult.objects.filter(Bin=None)
    orphan_verdict = PicaResult.objects.filter(Bin=None)
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
            clean_up_on_pipeline_fail(key, ppath, failtype="UNKNOWN")
            raise RuntimeError("Pipeline run has failed. Error status in DB has been updated.")

        remove_temp_files(infolder)

        # if no output file has been generated despite no pipeline error
        # (= if error checking consumes all input files)
        if not os.path.exists(os.path.join(outfolder, "{k}.zip".format(k=key))):
            clean_up_on_pipeline_fail(key, ppath, failtype="ALL_DROPPED")
            raise RuntimeError("No input files have passed error checking.")

        return pipeline_call.returncode
