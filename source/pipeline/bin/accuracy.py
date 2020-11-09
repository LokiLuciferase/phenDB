#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import Taxon, Bin, PicaModelAccuracy, PicaModel
import math

def round_nearest(x, a):
    return round(round(x / a) * a, -int(math.floor(math.log10(a))))


complecontaitem = sys.argv[1]
mdsum = sys.argv[2]
modelname = sys.argv[3]

# get completeness and contamination
try:
    with open(complecontaitem, "r") as ccfile:
        comple, conta, strainhet, taxid, tname, trank = ccfile.readline().strip().split("\t")
except ValueError:
    raise ValueError
    #comple, conta, strainhet, taxid, tname, trank = [0 for x in range(6)]

cc = [float(x) for x in [comple, conta, strainhet]]
for i in range(len(cc)):
    if cc[i] < 0:
        cc[i] = 0
    if cc[i] > 1:
        cc[i] = 1

# get taxonomic information
try:
    taxon_entry = Taxon.objects.get(tax_id=taxid)
    tid = taxon_entry.tax_id
except ObjectDoesNotExist:
    tid = "1"

# update bin with taxonomic and compleconta information
try:
    parentbin = Bin.objects.get(md5sum=mdsum)
    parentbin.comple = float(cc[0])
    parentbin.conta = float(cc[1])
    parentbin.strainhet = float(cc[2])
    parentbin.tax_id = tid
    parentbin.save()

except ObjectDoesNotExist:
    sys.exit("Bin not found.")

# check the balanced accuracy. other metrices would be very similar to evaluate, just change the part after the last dot
# the print statments includes a newline at the end, this is important for processing further downstream

print(PicaModelAccuracy.objects.get(model=PicaModel.objects.filter(model_name=modelname).latest('model_train_date'),
comple=round_nearest(float(cc[0]),0.05),
conta=round_nearest(float(cc[1]),0.05)).mean_balanced_accuracy)