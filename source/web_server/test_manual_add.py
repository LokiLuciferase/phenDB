#
# Created by Lukas Lüftinger on 08/11/2017.
#

# this script manually enters lines into database using django
# enters a job and an assorted bin into the database
# and prints them afterwards
#
import os
import django
from django.utils import timezone

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()

from phenotypePredictionApp.models import *

# newjob = job(user_ip="192.168.0.1", user_email="test@gmail.com", folder_path="/path/to/folder", job_status="FUBAR")
# newjob.save()
#
# newbin = bin(bin_id="abc",
#              file_name="abc.fasta",
#              job=newjob,
#              genome_path=r"D:\atav1st\Dropbox\coding\py\GitHub\phenDB\docs\TODO.txt",
#              md5sum="iamarandomstringlolol121212",
#              errors="FUUUUBAR",
#              comple=0.1,
#              conta=0.99999)
# newbin.save()
#
# newmodel = model(model_id="COOL", model_desc="Is the organism cool?", model_train_date=timezone.now())
# newmodel.save()
#
# newenog = enog(enog_id="SUNGLS", enog_descr="Sunglasses")
# newenog.save()
#
# newrank = model_enog_ranks(model=newmodel, enog=newenog, internal_rank=1)
# newrank.save()
#
# newresult_enog = result_enog(enog=newenog, bin=newbin)
# newresult_enog.save()
#
# newresult_model = result_model(bin=newbin, model=newmodel, verdict=True, accuracy=0.99)
# newresult_model.save()

## add job manually for testing purposes
samplejob = job(job_id="test_working", user_ip="192.168.0.1", user_email="abc@test.com", folder_path="/scratch/swe_ws17/data/test_working", job_status="should_work")
samplejob.save()


# print(job.objects.all())
# print(bin.objects.all())
# print(model.objects.all())
# print(enog.objects.all())
# print(model_enog_ranks.objects.all())
# print(result_model.objects.all())
# print(result_enog.objects.all())
