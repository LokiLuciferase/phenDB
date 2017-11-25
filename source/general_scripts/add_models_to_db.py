#!/usr/bin/env python3
import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *
import sys
from django.core.exceptions import ObjectDoesNotExist

def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1

def check_groupfile(enogname, enogrank):
    enog_rank_list_featurgroup=[]
    for groupline in groupfile:
        if enogname == groupline.split("\t")[0]:
            groupline = groupline.split("\t")[1]
            groupline = groupline.split("/")
            nr_of_enogs_in_fg= len(groupline)
            for entry in groupline:
                entry = entry.rstrip()
                try:
                    new_enog_rank = model_enog_ranks(model=newmodel,
                                                     enog=enog.objects.get(enog_name=entry),
                                                     internal_rank=enogrank/nr_of_enogs_in_fg)
                    enog_rank_list_featurgroup.append(new_enog_rank)
                    print(new_enog_rank)
                    # sys.stdout.write("Added Enog(s)+rank(s) contained in "
                    #                  " {gr}.\r".format(gr=enogname))
                    # sys.stdout.flush()
                except ObjectDoesNotExist:
                    print("\n cant find this ENOG:\n", entry)
            # if enog was a feature group contained in the rank.groups file
            #print()
            return enog_rank_list_featurgroup
    # if enog was something else
    print("\n cant find this ENOG or featuregroup:\n", enogname)
    return []

def rankfile_to_list():
    # make a list of enog_ranks which should be added to the db
    enog_rank_list=[]
    counter=0
    for line in rankfile.readlines()[1:]:
        line = line.split()
        counter+=1
        try:
            new_enog_rank = model_enog_ranks(model=newmodel, enog=enog.objects.get(enog_name=line[0]),
                                             internal_rank=float(line[1]))
            enog_rank_list.append(new_enog_rank)
            print(new_enog_rank)
            # sys.stdout.write("Added Enog+rank {el}.\r".format(el=line[0]))
            # sys.stdout.flush()
            sys.stdout.write("Added Enog Nr. {nr}     of {totnr}.\r".format(nr=counter, totnr=num_lines_ranksfile))
            sys.stdout.flush()

        except ObjectDoesNotExist:
            # if the enog does not exist in the annotationsdfile, it might be a "feature group"
            # check this by looking up in the .rank.groups file. If it is the case, add all enogs in the
            # feature group with the rank of the feature group to the db
            print("FG!!!!!!!!!!!!!!!!!!!!!!")
            add_this=check_groupfile(line[0], float(line[1]))
            if add_this:
                enog_rank_list+=add_this

    return enog_rank_list

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
    num_lines_ranksfile=file_len(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank")
    with open(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank","r") as rankfile:
        with open(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank.groups", "r") as groupfile:
            print("Creating list of enogs...")
            enog_rank_list_filled=rankfile_to_list()
            print(" \n Saving to database... ")
            model_enog_ranks.objects.bulk_create(enog_rank_list_filled)

print(model.objects.all())
