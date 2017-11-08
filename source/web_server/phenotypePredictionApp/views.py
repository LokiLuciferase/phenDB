from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from .forms import FileForm
import uuid
#from businessLogic.startProcess import demoProcess
from pprint import pprint

# Create your views here.

def index(request):
    print("index called")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResult' : 'none'}
    return HttpResponse(template.render(context, request))

def sendinput(request):
    print("sendinput called")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'There are your results',
               'showResult' : 'block'}
    postobj = request.POST.copy()
    fileobj = request.FILES.copy()

    #works only if just one file is uploaded
    for file_name, file in request.FILES.items():
        name = request.FILES[file_name].name
        postobj['file_name'] = name  #LL: indented this because I'm like 99% sure that's what you wanted

    form = FileForm(postobj, fileobj)
    if(form.is_valid()):
        print("is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        # startProcess(key)
    return HttpResponse(template.render(context, request))
