# !/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import Bin, Job, HmmerResult, BinInJob

jobname = sys.argv[1]
nr_of_files = int(sys.argv[2])
mdsum = sys.argv[3]
binname = sys.argv[4]

try:
    parentjob = Job.objects.get(key=jobname)
    parentjob.total_bins = nr_of_files
    parentjob.save()

except ObjectDoesNotExist:
    sys.exit("Job not found.")

# if the bin exists, calculate its hmmer and compleconta file from the db entries
try:
    thisbin = Bin.objects.get(md5sum=mdsum)

    with open("reconstructed_hmmer_file.txt", "w") as hmmer:
        for entry in HmmerResult.objects.filter(bin=thisbin):
            hmmer.write("dummy\t" + entry.enog.enog_name + "\tdummy\n")

    with open("reconstructed_compleconta_file.txt", "w") as complecon:
        recon_cc_string = "\t".join([str(x) for x in (thisbin.comple,
                                                       thisbin.conta,
                                                       thisbin.strainhet,
                                                       thisbin.tax_id,
                                                       "dummy", "dummy")])
        complecon.write(recon_cc_string)

    print("NO", end='')

    # Also update the job completeness here, because then hmmer wont be called
    current_finished = parentjob.finished_bins
    total_jobs = parentjob.total_bins
    if current_finished < total_jobs - 1:
        parentjob.finished_bins = int(current_finished) + 1
        parentjob.save()


except ObjectDoesNotExist:
    # write impossible Nr. for comple and conta that would cause an error during get_accuracy if not overwritten:
    thisbin = Bin(md5sum=mdsum, comple=2, conta=2, strainhet=2)
    thisbin.save()
    print("YES", end='')

    # if it not in the db yet, just create dummy files
    with open("reconstructed_hmmer_file.txt", "w") as hmmer:
        hmmer.write("dummy")
    with open("reconstructed_compleconta_file.txt", "w") as complecon:
        complecon.write("dummy")

try:
    assoc = BinInJob(bin=thisbin, job=parentjob, bin_alias=binname)
    assoc.save()
except IntegrityError:
    if binname.startswith("PHENDB_PRECALC"):  # allow muliple identical files in same job if precalculation job
        pass
    else:
        raise