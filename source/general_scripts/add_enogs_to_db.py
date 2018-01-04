import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import gzip
import sys

ENOG_files = ["/var/www/eggnog/eggnog_4.5/data/bacNOG/bacNOG.annotations.tsv.gz",
              "/var/www/eggnog/eggnog_4.5/data/bactNOG/bactNOG.annotations.tsv.gz" ]

enoglist=[]
namlist=[]
for annot in ENOG_files:
    print("Processing ", annot)
#Open the annotations file, every line is an enog. Translate the Lettercode in the annotations-file using the
#descriptions file and save the enog + description + categories to the db.
    with gzip.open(annot,mode="rt") as f:
        counter=0
        for line in f:
            line=line.split("\t")
            enoglist.append(enog(enog_name=line[1], enog_descr=line[5].rstrip() ))
            namelist.append(line[1])
            counter+=1
            sys.stdout.write('\r')
            sys.stdout.write("Adding enog nr. ")
            sys.stdout.write(str(counter))
            sys.stdout.flush()

print("\n writing enogs to db...")

print(len(namelist) != len(set(namelist)))
enog.objects.bulk_create(enoglist)