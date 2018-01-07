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
from ipware.ip import get_real_ip

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
    postobj['user_ip'] = get_real_ip(request)
    form = FileForm(postobj, fileobj)
    if(form.is_valid()):
        print("Is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        thread = startProcessThread(key)
        thread.start()
    resultObj = UploadedFile.objects.get(key=key)
    return redirect(resultObj)

def getResults(request):
    #TODO: show error on UI not related to sanity-error
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    key = getKeyFromUrl(request)
    obj = UploadedFile.objects.get(key=key)
    #errors: when errors from model UploadedFile set to true -> sanity error, when None -> Error in the pipeline
    error_severity = 'warn'
    if obj.job_status == '100':
        showResult = 'block'
        showProgressBar = False
        refresh = False
        if(obj.errors == None):
            showErrorMessage = True
            error_severity = 'error'
            error_summary = 'unknown error'
            error_message = 'phendb service not working'
        else:
            showErrorMessage = obj.errors
            error_summary = 'Invalid input file(s)'
            error_message = 'Please check the invalid_input_files.log file in the summaries directory'
    else:
        showResult = 'none'
        showProgressBar = True
        refresh = True
        showErrorMessage = False
        error_message = ""
        error_summary = ""
        error_severity = ""
    print('status ' + obj.job_status)
    context = {'result' : 'download/',
               'showResult' : showResult,
               'showProgressBar' : showProgressBar,
               'progress' : obj.job_status,
               'refresh' : refresh,
               'showInputForm': 'none',
               'showErrorMessage': showErrorMessage,
               'error_severity' : error_severity,
               'error_summary' : error_summary,
               'error_message' : error_message}
    return HttpResponse(template.render(context, request))

def fileDownload(request):
    key = getKeyFromUrl(request)
    resFile = UploadedFile.objects.get(key = key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="{k}.zip"'.format(k=key)
    return response