import time
import argparse
import os
import django
from django.core.files import File
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *

def main():
    print("fakeScript called")
    parser = argparse.ArgumentParser()
    parser.add_argument("--infolder")
    parser.add_argument("--key")
    args = parser.parse_args()
    time.sleep(60)
    obj = UploadedFile.objects.filter(key=args.key)
    file = open(args.infolder, 'rb')
    djangoFile = File(file)
    obj[0].fileOutput.save('result.tar.gz', djangoFile, save="True")
    obj.update(job_status = '100')


main()