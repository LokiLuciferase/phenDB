#!/usr/bin/env python3

import sys
import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import Bin, PicaModel, Enog, PicaResultExplanation

mdsum_file = sys.argv[1]
mdsum = sys.argv[2]

DB_ENOGS = {v.enog_name: v for k, v in Enog.objects.in_bulk().items()}


# get Bin from db
try:
    parentbin = Bin.objects.get(md5sum=mdsum)
except ObjectDoesNotExist:
    sys.exit("Bin not found.")

conditions = []
with open(mdsum_file, "r") as picaresults:
    for line in picaresults:
        conditions.append(line.split())

expl_list = []
for result in conditions:  # result = [modelname, rank, enog_name, is_present, delta_shap]
    try:
        model_name, explanation_rank, enog_name, enog_is_present, delta_shap = (
            str(result[0]),
            int(result[1]),
            str(result[2]),
            bool(float(result[3])),
            float(result[4]),
        )
        # get model from db
        try:
            this_model = PicaModel.objects.filter(model_name=model_name).latest('model_train_date')
        except ObjectDoesNotExist:
            sys.exit("Current Model for this result not found.")

        expl = PicaResultExplanation(
            model=this_model,
            bin=parentbin,
            enog=DB_ENOGS[enog_name],
            enog_is_present=enog_is_present,
            delta_shap=delta_shap
        )
        expl_list.append(expl)
    except (IntegrityError,) as e:
        sys.exit(e)

try:
    PicaResultExplanation.objects.bulk_create(expl_list)
# TODO: inconsistency of pica pval despite same model version??
except IntegrityError:
    for result in expl_list:
        if not PicaResultExplanation.objects.filter(
            model=result.model,
            bin=result.bin,
            enog=result.enog
        ).exists():
            result.save()
            print(result)
