import os
import django
from django.utils import timezone
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *



PICAMODELFOLDER="/Users/peterpeneder/Documents/UNI-BIOINFO/Softwareentwicklung/models2"
all_picamodels=os.listdir(PICAMODELFOLDER)

for picamodel in all_picamodels:
    #TODO: change model_train_date?

    if model.objects.filter(model_id=picamodel).exists():
        # caclulate the new version nr from the highest version nr of already existing entries for the model
        new_version_nr = 1 + max(existing_model.version_nr for existing_model in model.objects.filter(model_id=picamodel))
        # set the previously most current version to non-newest
        model.objects.filter(model_id=picamodel, is_newest=True).update(is_newest=False)
    else:
        new_version_nr=1
    newmodel = model(model_id=picamodel, model_desc="", is_newest=True, version_nr=new_version_nr, model_train_date=timezone.now())
    newmodel.save()

    #TODO: parallelize this?
    with open(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank","r") as rankfile:
        counter=0
        for line in rankfile:
            if counter!=0:
                line=line.split()
                new_enog = enog(enog_id=line[0], enog_descr="")
                new_enog_rank = model_enog_ranks(model=newmodel, enog=new_enog, internal_rank=line[1])
                new_enog_rank.save()
                new_enog.save()
            counter+=1

print(model.objects.all())
print(enog.objects.all())
print(model_enog_ranks.objects.all())

#model.objects.all().delete()

# newmodel = model(model_id=picamodel, model_desc="", version_nr=3, model_train_date=timezone.now())
# newmodel.save()