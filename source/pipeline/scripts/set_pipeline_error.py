#
# Created by Lukas Lüftinger on 04/01/2018.
#
import os
import sys
import django
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import UploadedFile

jobkey = sys.argv[1]

current_job = UploadedFile.objects.get(key=jobkey)
current_job.errors = None
current_job.save()
