#
# Created by Lukas Lüftinger on 6/29/18.
#
import os
import sys

import django

ppath = "/apps/phenDB/source/web_server"
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
sys.path.append(ppath)
django.setup()

from phenotypePredictionApp.models import Bin, Job, BinInJob, Taxon, PicaModel, PicaModelTrainingData

def read_files(path):
    fd = {}
    for file in os.listdir(path):
        modelname = file.replace(".training.txt", "")
        fullname = os.path.join(path, file)
        with open(fullname, "r") as infile:
            infile.readline()
            lines = [x.strip().split("\t") for x in infile.readlines()]
        lines_reordered = [(x[2], x[0], True if x[1] == "YES" else False) for x in lines]
        fd[modelname] = lines_reordered
    return fd

uld = read_files("/apps/PICA/training_files_accessions")

for mname in uld.keys():
    oldmodel = PicaModel.objects.filter(model_name=mname).latest("model_train_date")
    print(oldmodel)
    training_data = []
    for row in uld[mname]:
        print(row)
        training_data.append(PicaModelTrainingData(model=oldmodel,
                                                   tax_id=str(row[0]),
                                                   assembly_id=str(row[1]),
                                                   verdict=str(row[2])))
    try:
        PicaModelTrainingData.objects.bulk_create(training_data)
    except Exception as e:
        print(e)

