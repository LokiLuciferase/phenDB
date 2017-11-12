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
    newmodel = model(model_desc="picamodel", version_nr=3, model_train_date=timezone.now())
    newmodel.save()
 #   if model.objects.filter(picamodel in model_id).exists():
  #      print(model.objects.filter(model_id=picamodel))
        #print(model.objects.get(model_id=picamodel).version_nr)
        #new_version_nr=1+model.objects.filter(model_id=picamodel)
    #     else:
    #     new_version_nr=1
    # newmodel = model(model_id=picamodel, model_desc="", version_nr=new_version_nr, model_train_date=timezone.now())
    # newmodel.save()
    # with open(PICAMODELFOLDER+"/"+picamodel+"/"+picamodel+".rank","r") as rankfile:
    #     counter=0
    #     for line in rankfile:
    #         if counter==0:
    #             print (line)
    #         counter+=1

print(model.objects.all())
#model.objects.all().delete()