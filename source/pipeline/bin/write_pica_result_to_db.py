# !/usr/bin/env python3

import sys
import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import Job, Bin, PicaModel, PicaResult

mdsum_file = sys.argv[1]
mdsum = sys.argv[2]

# get Bin from db
try:
    parentbin = Bin.objects.get(md5sum=mdsum)
except ObjectDoesNotExist:
    sys.exit("Bin not found.")

conditions = []
with open(mdsum_file, "r") as picaresults:
    for line in picaresults:
        conditions.append(line.split())

modelresultlist = []
for result in conditions:  # result = [modelname, verdict, pica_p_val, balanced_accuracy]
    try:
        get_bool = {"YES": True, "NO": False, "N/A": None}
        if result[2] == "NA" or result[2] == "N/A":
            get_pica_pval = float(0)
        else:
            get_pica_pval = float(result[2])
        boolean_verdict = get_bool[result[1]]
        # get model from db
        try:
            this_model = PicaModel.objects.filter(model_name=result[0]).latest('model_train_date')
        except ObjectDoesNotExist:
            sys.exit("Current Model for this result not found.")

        modelresult = PicaResult(verdict=boolean_verdict,
                                 accuracy=float(result[3]),
                                 pica_pval=get_pica_pval,
                                 bin=parentbin,
                                 model=this_model
                                 )
        modelresultlist.append(modelresult)
    except (IntegrityError,) as e:
        sys.exit(e)
try:
    PicaResult.objects.bulk_create(modelresultlist)
except IntegrityError:
    for result in modelresultlist:
        PicaResult.objects.get_or_create(verdict=result.verdict,
                                         accuracy=result.accuracy,
                                         pica_pval=result.pica_pval,
                                         bin=result.bin,
                                         model=result.model)
