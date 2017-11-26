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
    if(reqPath[-1] == '/'):
        reqPath = reqPath[:-1]
    urlparts = reqPath.rsplit('/')
    for part in urlparts:
        if len(part) == 36:
            return part
    return None



#-------------------Views-------------------------------------------------

def index(request):
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResult' : 'none',
               'showInputForm': 'block',
               'showProgressBar' : False,
               'refresh' : False,}
    return HttpResponse(template.render(context, request))

def sendinput(request):
    postobj = request.POST.copy()
    fileobj = request.FILES.copy()
    key = str(uuid.uuid4())
    postobj['key'] = key

    #works only if just one file is uploaded
    for filename, file in request.FILES.items():
        name = request.FILES[filename].name
    postobj['filename'] = name
    postobj['fileInput'] = fileobj['fileInput']
    form = FileForm(postobj, fileobj)
    if(form.is_valid()):
        print("is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        thread = startProcessThread(key)
        thread.start()
    resultObj = UploadedFile.objects.get(key=key)
    return redirect(resultObj)

def getResults(request):
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    key = getKeyFromUrl(request)
    obj = UploadedFile.objects.get(key=key)
    if obj.job_status == '100':
        showResult = 'block'
        showProgressBar = False
        refresh = False
    else:
        showResult = 'none'
        showProgressBar = True
        refresh = True
    print('status ' + obj.job_status)
    #TODO: !!!! ERROR MESSAGE NOT WORKING !!!
    context = {'result' : 'download/',
               'showResult' : showResult,
               'showProgressBar' : showProgressBar,
               'progress' : obj.job_status,
               'refresh' : refresh,
               'showInputForm': 'none',
               'showErrorMessage': True,
               'error_message' : 'Test'}
    return HttpResponse(template.render(context, request))

def fileDownload(request):
    key = getKeyFromUrl(request)
    resFile = UploadedFile.objects.get(key = key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="phendb_results.tar.gz"'
    return response