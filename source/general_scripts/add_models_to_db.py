#!/usr/bin/env python3
import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import sys
from django.core.exceptions import ObjectDoesNotExist


#PICAMODELFOLDER="/Users/peterpeneder/Desktop/models/models2"
PICAMODELFOLDER="/scratch/swe_ws17/data/models"
all_picamodels=os.listdir(PICAMODELFOLDER)
for picamodel in all_picamodels:

    #TODO: change model_train_date?
    #check if there alreads is a model with the same name in the db, if yes, add new entry with higher version
    #nr and set all others to non-active
    if model.objects.filter(model_name=picamodel).exists():
        sys.stdout.write("Updating model {m}.\n".format(m=picamodel))
        sys.stdout.flush()

        # caclulate the new version nr from the highest version nr of already existing entries for the model
        new_version_nr = 1 + max(existing_model.version_nr for existing_model in model.objects.filter(model_name=picamodel))
        # set the previously most current version to non-newest
        model.objects.filter(model_name=picamodel, is_newest=True).update(is_newest=False)

    else:
        sys.stdout.write("#########\nAdding new model {pm}.\n########\n".format(pm=picamodel))
        sys.stdout.flush()
        new_version_nr=1

    try: #check if there is a .description file, if there is, read the description
        with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".description", "r") as descfile:
            desc=""
            for line in descfile:
                desc+=" "+line.rstrip()
    except FileNotFoundError:
        desc=""

    # add the model to the database
    newmodel = model(model_name=picamodel, model_desc=desc, is_newest=True, version_nr=new_version_nr,
                     model_train_date=timezone.now())
    newmodel.save()



    # add the ranks of the models to the database
    #TODO: parallelize this?
    #read the .rank file of the model and extract enogs and their ranks
    with open(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank","r") as rankfile:
        with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".rank.groups", "r") as groupfile:
            grouplist = groupfile.readlines()
            for line in rankfile.readlines()[1:]:
                #try:
                line=line.split()
                try:
                    new_enog_rank = model_enog_ranks(model=newmodel, enog=enog.objects.get(enog_name=line[0]),
                                                     internal_rank=line[1])
                    new_enog_rank.save()
                    sys.stdout.write("Added Enog+rank {el}.\r".format(el=line[0]))
                    sys.stdout.flush()
                except ObjectDoesNotExist:
                    #if the enog does not exist in the annotationsdfile, it might be a "feature group"
                    # check this by looking up in the .rank.groups file. If it is the case, add all enogs in the
                    # feature group with the rank of the feature group to the db

                    for groupline in grouplist:
                        if line[0] == groupline.split("\t")[0]:
                            groupline=groupline.split("\t")[1].split("/")
                            for entry in groupline:
                                entry = entry.rstrip()
                                try:
                                    new_enog_rank = model_enog_ranks(model=newmodel,
                                                                     enog=enog.objects.get(enog_name=entry),
                                                                     internal_rank=line[1])
                                    new_enog_rank.save()
                                    sys.stdout.write("Added Enog(s)+rank(s) contained in "
                                                     " {gr}.\r".format(gr=line[0]))
                                    sys.stdout.flush()
                                except ObjectDoesNotExist:
                                    print("\n cant find this ENOG:\n", entry)
                        # if enog was a feature group contained in the rank.groups file
                        print()
                    # if enog was something else
                    print("\n cant find this ENOG or featuregroup:\n", line[0])


                    # todo: exception for when its not a feature group

                                        #except:
                 #   sys.stdout.write("Skipping.\n")
                  #  pass

print(model.objects.all())
