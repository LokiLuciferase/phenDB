#!/usr/bin/env python3

import os
import django
import json
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import sys
from django.core.exceptions import ObjectDoesNotExist


PICAMODELFOLDER="/scratch/swe_ws17/data/models"
all_picamodels=os.listdir(PICAMODELFOLDER)
for picamodel in all_picamodels:
    jsonfile=PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".accuracy.json"
    data = json.load(open(jsonfile))
    acc_list=[]
    for i in range (0,21*21):
        this_data=data[i]
        acc_list.append(model_accuracies(model=model.objects.get(model_name=picamodel, is_newest=True),
                                       comple=this_data["completeness"], conta=this_data["contamination"],
                                       mean_balanced_accuracy=this_data["mean_balanced_accuracy"],
                                       mean_fp_rate = this_data["mean_fp_rate"],
                                       mean_fn_rate = this_data["mean_fn_rate"]))
    # for comple in range(0,21):
    #     for conta in range(0, 21):
    #         print("given: ", comple*0.05, conta*0.05)
    #         print("looked comple", data[conta+comple*21]["completeness"], " looked conta", data[conta+comple*21]["contamination"])
    #         print()
    #

