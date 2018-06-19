#!/usr/bin/env python3
import os
import sys
import timeit
import json

import django
from django.core.exceptions import ObjectDoesNotExist
from django.utils import timezone

django.setup()
from phenotypePredictionApp.models import *


PICAMODELFOLDER = "/apps/PICA/models"
all_picamodels = os.listdir(PICAMODELFOLDER)


def file_len(fname):
    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def rankfile_to_list(rankfile, groupfile, db_enogs):
    # make a list of enog_ranks which should be added to the db
    enog_rank_list = []

    # create a dictionary over all entries in the .ranks.group file, with the featuregroup-name as key and as list of
    # corresponding enogs as value
    featuregroup_dict = {}
    for groupline in groupfile:
        groupname, enogs_in_group = groupline.rstrip().split("\t")
        enogs_in_group = enogs_in_group.split("/")
        featuregroup_dict[groupname] = enogs_in_group

    # skip first line
    for counter, line in enumerate(rankfile):
        if counter == 0:
            continue
        enog_name, score, verdict = line.split()

        try:
            # check if the enog is contained in the database, if yes, add the enog_ranks object to the enog_rank_list try:
            new_enog_rank = EnogRank(model=newmodel,
                                     enog=db_enogs[enog_name],
                                     internal_rank=counter,
                                     score=float(score),
                                     pred_class=(False if float(score) <= 0 else True))
            enog_rank_list.append(new_enog_rank)
            sys.stdout.write("Added Enog Nr. {nr}        of {totnr}.\r".format(nr=counter, totnr=num_lines_ranksfile))
            sys.stdout.flush()
        # If the "enog" is not contained in the db, it might actually be a "feature group"
        # check this by looking up in the feature_group dict. If it is the case, add all enogs in the
        # feature group to the db, each with a rank of "rank of fg"/"number of enogs in fg"
        except KeyError:
            try:
                nr_of_enogs_in_fg = len(featuregroup_dict[enog_name])
                for fg_enog in featuregroup_dict[enog_name]:
                    try:
                        new_enog_rank = EnogRank(model=newmodel,
                                                 enog=db_enogs[fg_enog],
                                                 internal_rank=counter,
                                                 score=float(score) / nr_of_enogs_in_fg,
                                                 pred_class=(False if float(score) <= 0 else True))
                        enog_rank_list.append(new_enog_rank)
                    except KeyError:
                        sys.exit("\n ERROR: The .rank.groups file of the model contained {fg_enog} which could not "
                                 "be found in the database. ABORTING. \n".format(fg_enog=fg_enog))

            except KeyError:
                sys.exit("\n ERROR. The .rank file of the model contained {enog_name} which could be found neither "
                         "in the database nor in the .rank.groups file. ABORTING. \n".format(enog_name=enog_name))
    return enog_rank_list

# store a dictionary of all enogs that are currently in the db
# if statement forces the objects to be loaded from the database (?)
db_enogs = Enog.objects.in_bulk()
db_enogs = {v.enog_name: v for k, v in db_enogs.items()}

if db_enogs:

    for picamodel in all_picamodels:
        # check if there alreads is a model with the same name in the db, if yes, add new entry with higher version
        # nr and set all others to non-active
        if PicaModel.objects.filter(model_name=picamodel).exists():
            sys.stdout.write("#########\nUpdating model {m}:\n\n".format(m=picamodel))
            sys.stdout.flush()

            # caclulate the new version nr from the highest version nr of already existing entries for the model
        # new_version_nr = 1 + max(existing_model.version_nr for existing_model in model.objects.filter(model_name=picamodel))
        # set the previously most current version to non-newest
        # model.objects.filter(model_name=picamodel, is_newest=True).update(is_newest=False)

        else:
            sys.stdout.write("\n#########\nAdding new model {pm}:\n\n".format(pm=picamodel))
            sys.stdout.flush()
        # new_version_nr=1

        try:  # check if there is a .description file, if there is, read the description
            with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".description", "r") as descfile:
                desc = ""
                for line in descfile:
                    desc += " " + line.rstrip()
        except FileNotFoundError:
            desc = ""

        try:  # check if there is a .type file, if there is, read the description
            with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".type", "r") as typefile:
                for line in typefile:
                    type += line.rstrip()
        except FileNotFoundError:
            type = "NA"

        # add the model to the database
        newmodel = PicaModel(model_name=picamodel, model_desc=desc, type=type,
                             model_train_date=timezone.now())
        newmodel.save()

        # add the ranks of the models to the database
        # read the .rank file of the model and extract enogs and their ranks
        num_lines_ranksfile = file_len(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".rank")
        with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".rank", "r") as rankfile:
            with open(PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".rank.groups", "r") as groupfile:
                print("Creating list of enogs...")
                try:
                    enog_rank_list_filled = rankfile_to_list(rankfile, groupfile, db_enogs)
                except:  # specifc error here?
                    newmodel.delete()
                    sys.exit("\n There was a problem when creating enog lists")

                print(" \n Saving to database... ")
                try:
                    EnogRank.objects.bulk_create(enog_rank_list_filled)
                except Exception as e:  # specifc error here?
                    newmodel.delete()
                    print(e)
                    sys.exit("\n There was a problem when writing enog_ranks-objects to the db")

        print("Processing the model's accuracy file...")
        # Read the model's accuracy file and enter the data into the db
        try:
            jsonfile = PICAMODELFOLDER + "/" + picamodel + "/" + picamodel + ".accuracy.json"
            data = json.load(open(jsonfile))
            acc_list = []
            for i in range(0, 21 * 21):
                this_data = data[i]
                acc_list.append(PicaModelAccuracy(model=newmodel,
                                                  comple=this_data["completeness"], conta=this_data["contamination"],
                                                  mean_balanced_accuracy=this_data["mean_balanced_accuracy"],
                                                  mean_fp_rate=this_data["mean_fp_rate"],
                                                  mean_fn_rate=this_data["mean_fn_rate"]))
            print("writing accuracy entries into db...")
            PicaModelAccuracy.objects.bulk_create(acc_list)

        except Exception as e:  # specifc error here?
            newmodel.delete()
            print(e)
            sys.exit("\n There was a problem while processing the accuracy file or while writing it to the db ")

        print("Attempting to add information on training data of model...")
        # check if the model has a .phenotype file, and upload information to model_used_genomes

        try:
            train_list = []
            phenotypefile = os.path.join(PICAMODELFOLDER, picamodel, picamodel + ".phenotype")
            taxidfile = os.path.join(PICAMODELFOLDER, picamodel, picamodel + ".taxids")
            with open(phenotypefile, "r") as pin:
                accession_type, trait = pin.readline().strip().split("\t")
                train_data = [x.strip().split() for x in pin.readlines()]
            with open(taxidfile, "r") as tin:
                taxid_list = [x.strip().split("\t") for x in tin.readlines()]
                taxid_dict = {x: y for x, y in taxid_list}

            train_list = [PicaModelTrainingData(model=newmodel,
                                                tax_id=taxid_dict[x],
                                                assembly_id=x,
                                                verdict=(True if y == "YES" else False)) for x, y in train_data]
            print("Writing training data to db...")
            PicaModelTrainingData.objects.bulk_create(train_list)

        except FileNotFoundError:
            print("No training data was found for this model. Skipping.")
