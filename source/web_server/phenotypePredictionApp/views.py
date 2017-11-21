from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .forms import FileForm
from django.shortcuts import redirect
import uuid
from businessLogic.startProcess import startProcess
from phenotypePredictionApp.models import UploadedFile
from pprint import pprint

# Create your views here.

def index(request):
    print("index called")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResult' : 'none',
               'showProgressBar' : 'none',
               'refresh' : False}
    return HttpResponse(template.render(context, request))

def sendinput(request):
    print("sendinput called")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'There are your results',
               'showResult' : 'block'}
    postobj = request.POST.copy()
    fileobj = request.FILES.copy()
    key = str(uuid.uuid4())
    postobj['job_name'] = key

    #works only if just one file is uploaded
    for filename, file in request.FILES.items():
        name = request.FILES[filename].name
    postobj['filename'] = name
    postobj['upload_path'] = fileobj['upload_path']

    print(postobj['filename'])
    print(postobj['job_name'])


    form = FileForm(postobj, fileobj)

    print(form.data)
    print('\n')
    print(fileobj)

    if(form.is_valid()):
        print("is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        startProcess(key)

    resultObj = UploadedFile.objects.get(key=key)
    return redirect(resultObj)

def getResults(request):
    print("works")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResult' : 'none',
               'showProgressBar' : 'block',
               'refresh' : True}
    return HttpResponse(template.render(context, request))