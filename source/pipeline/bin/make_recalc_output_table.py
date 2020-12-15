#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 10/2/18.
#
import sys
import os
import math

import django
import numpy as np
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError


def round_nearest(x, a):
    return round(round(x / a) * a, -int(math.floor(math.log10(a))))


os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import (
    Bin,
    BinInJob,
    Job,
    PicaModel,
    HmmerResult,
    Enog,
    PicaModelAccuracy,
    PicaResult,
)

return_tuples = []
precalc_bins = list(Bin.objects.filter(bininjob__job__key__icontains="PHENDB_PRECALC"))
precalc_bins = sorted(precalc_bins, key=lambda x: x.md5sum)
model_names = [x["model_name"] for x in list(PicaModel.objects.values("model_name").distinct())]
batch_number = int(sys.argv[2])
total_batch_no = int(sys.argv[3])
if batch_number >= 0:
    split_bins = np.array_split(precalc_bins, total_batch_no)
    precalc_bins = list(split_bins[batch_number])

for bin in precalc_bins:
    hmmerpath = "{md5sum}_RECALC.out".format(md5sum=bin.md5sum)
    enoglist = [
        x["enog_name"] for x in Enog.objects.filter(hmmerresult__bin=bin).values("enog_name")
    ]
    write_out = "\n".join(["\t".join(["dummy", x, "dummy"]) for x in enoglist])
    with open(hmmerpath, "w") as hmmerfile:
        hmmerfile.write(write_out)
    for model_name in model_names:
        newest_model = PicaModel.objects.filter(model_name=model_name).latest("model_train_date")
        if PicaResult.objects.filter(
            bin=bin, model=newest_model
        ).exists():  # don't recalc if for same bin and same model results exist
            continue
        binname = bin.md5sum
        md5sum = bin.md5sum
        model = "/".join([sys.argv[1], model_name])
        hmmeritem = os.path.abspath(hmmerpath)
        accuracy = PicaModelAccuracy.objects.get(
            model=newest_model,
            comple=round_nearest(float(bin.comple), 0.05),
            conta=round_nearest(float(bin.conta), 0.05),
        ).mean_balanced_accuracy
        return_tuples.append((binname, md5sum, model, hmmeritem, str(accuracy)))

if not return_tuples:
    return_tuples.append(("dummy_bin", "dummy_md5sum", "dummy_model", "dummy_hmmer", "dummy_acc"))

with open("recalc_table.csv", "w") as out_table:
    for t in return_tuples:
        out_table.write("\t".join(t))
        out_table.write("\n")
