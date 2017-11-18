import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import gzip

#path to the NOG annotations file: location on server:
#annot="/mirror/eggnog/eggnog_4.5/data/NOG/NOG.annotations.tsv.gz"
#for local testing:
annot="/Users/peterpeneder/Documents/UNI-BIOINFO/Softwareentwicklung/NOG.annotations.tsv.gz"

with gzip.open(annot,mode="rt") as f:
    for line in f:
        line=line.split("\t")
        catlist=[]
        with open ("/Users/peterpeneder/Documents/UNI-BIOINFO/Softwareentwicklung/NOG_descriptions.txt", 'r') as cat:
            for my_letter in line[4]:
                for catline in cat:
                    if "["+my_letter+"]" in catline:
                        catlist.append(catline)
        #if len(catlist)>1:
        print(catlist)
        new_enog = enog(enog_id=line[1], enog_descr=line[5])


