from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .forms import FileForm
from django.shortcuts import redirect
import uuid
from businessLogic.startProcess import startProcess
from phenotypePredictionApp.models import job
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
    job_name = str(uuid.uuid4())
    postobj['job_name'] = job_name

    #works only if just one file is uploaded
    for filename, file in request.FILES.items():
        name = request.FILES[filename].name
    postobj['filename'] = name
    postobj['upload_path'] = fileobj['fileInput']

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

        newjob = job(job_name=job_name,
                     filename=postobj['filename'],
                     )
        newjob.save()
    else:
        print("Form is not valid!")

        startProcess(job_name)

    resultObj = job.objects.get(job_name=job_name)
    return redirect(resultObj)

def getResults(request):
    print("works")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResult' : 'none',
               'showProgressBar' : 'block',
               'refresh' : True}
    return HttpResponse(template.render(context, request))