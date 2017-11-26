import time
import argparse
import os
import django
from django.core.files import File
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--infolder")
    parser.add_argument("--key")
    args = parser.parse_args()
    obj = UploadedFile.objects.filter(key=args.key)
    x = 0
    while x<6:
        time.sleep(1)
        val = 10*x
        print('val ' + str(val))
        obj.update(job_status=str(val))
        x += 1

    file = open(args.infolder, 'rb')
    djangoFile = File(file)
    obj[0].fileOutput.save('result.tar.gz', djangoFile, save="True")
    time.sleep(1)
    obj.update(job_status = '100')



main()

