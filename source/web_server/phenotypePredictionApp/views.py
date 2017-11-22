from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django.urls import reverse
from .forms import FileForm
from django.shortcuts import redirect
import uuid
from django.core.urlresolvers import resolve
from businessLogic.startProcess import startProcessThread
from phenotypePredictionApp.models import UploadedFile
from pprint import pprint

#------------------functions---------------------------------------------
#useful functions, NO views

def getKeyFromUrl(request):
    reqPath = request.path
    print(reqPath)
    if(reqPath[-1] == '/'):
        reqPath = reqPath[:-1]
    urlparts = reqPath.rsplit('/')
    for part in urlparts:
        if len(part) == 36:
            return part
    return None



#-------------------Views-------------------------------------------------

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
    postobj['key'] = key

    #works only if just one file is uploaded
    for filename, file in request.FILES.items():
        name = request.FILES[filename].name
    postobj['filename'] = name
    postobj['fileInput'] = fileobj['fileInput']

    print(postobj['filename'])
    print(postobj['key'])


    form = FileForm(postobj, fileobj)

    print(form.data)
    print('\n')
    print(fileobj)

    if(form.is_valid()):
        print("is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        thread = startProcessThread(key)
        thread.start()

    #TESTING -> DO NOT USE
    #res = ResultFile.objects.create(actualID=key, document=fileobj['fileInput'])

    resultObj = UploadedFile.objects.get(key=key)
    return redirect(resultObj)

def getResults(request):
    print("works")
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    key = getKeyFromUrl(request)
    obj = UploadedFile.objects.get(key=key)
    if obj.job_status == '100':
        showResult = 'block'
        refresh = 'false'
    else:
        showResult = 'none'
        refresh = 'true'
    context = {'result' : 'download/',
               'showResult' : showResult,
               'showProgressBar' : 'block',
               'progress' : obj.job_status,
               'refresh' : refresh}
    return HttpResponse(template.render(context, request))

def fileDownload(request):
    print("fileDownload called")
    key = getKeyFromUrl(request)
    resFile = UploadedFile.objects.get(key = key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="phendb_results.tar.gz"'
    return response