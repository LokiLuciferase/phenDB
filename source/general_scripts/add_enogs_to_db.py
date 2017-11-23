import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import gzip
import sys

# this does not work in VM! no eggnog mirror in vm!
ENOG_files = ["/mirror/eggnog/eggnog_4.5/data/bacNOG/bacNOG.annotations.tsv.gz",
              "/mirror/eggnog/eggnog_4.5/data/bactNOG/bactNOG.annotations.tsv.gz" ]

for annot in ENOG_files:
    print("Processing ", annot)
#Open the annotations file, every line is an enog. Translate the Lettercode in the annotations-file using the
#descriptions file and save the enog + description + categories to the db.
    with gzip.open(annot,mode="rt") as f:
        counter=0
        for line in f:
            line=line.split("\t")
            catlist=[]
            with open ("/mirror/eggnog/eggnog_4.5/COG_functional_categories.txt", 'r') as cat:
                for catline in cat:
                    for my_letter in line[4]:
                        if "["+my_letter+"]" in catline:
                            catlist.append(catline.rstrip())
            catlist = "; ".join(catlist)
            new_enog = enog(enog_name=line[1], enog_descr=line[5].rstrip(), enog_category=catlist)
            new_enog.save()
            counter+=1
            sys.stdout.write('\r')
            sys.stdout.write("Adding enog nr. ")
            sys.stdout.write(str(counter))
            sys.stdout.flush()

