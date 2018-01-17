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
    obj = UploadedFile.objects.filter(key=args.key)
    '''file = open(args.infolder, 'rb')
    djangoFile = File(file)
    obj[0].fileOutput.save('result.tar.gz', djangoFile, save="True") '''
    obj.update(total_bins = 5)
    counter = 0
    while(counter < 5):
        time.sleep(10)
        counter += 1
        obj.update(finished_bins = counter)


main()