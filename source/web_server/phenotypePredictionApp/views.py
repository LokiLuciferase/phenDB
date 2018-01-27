from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django.urls import reverse
from .forms import FileForm
from django.shortcuts import redirect
from businessLogic.mailNotification import MailNotification
import uuid
from django.core.urlresolvers import resolve
from businessLogic.startProcess import StartProcessThread
from phenotypePredictionApp.models import UploadedFile
from pprint import pprint
from ipware.ip import get_real_ip
from redis import Redis
from rq import Queue, get_current_job

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

def get_current_position(keyname):
    try:
        que = Queue('phenDB', connection=Redis())
        jobs = list(que.jobs)
        for index, job in enumerate(jobs):
            if job.id == keyname:
                return (index, len(jobs))
        return (-1, len(jobs))
    except:
        return (None, None)


#-------------------Views-------------------------------------------------

def index(request):
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResultCSS' : 'none',
               'showNotification' : False,
               'showInputFormCSS': 'block',
               'showProgressBar' : False,
               'refresh' : False}
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
    postobj['errors'] = False
    form = FileForm(postobj, fileobj)
    if(form.is_valid()):
        print("Is valid")
        modelInstance = form.save(commit=False)
        modelInstance.save()
        StartProcessThread(key).start()
        if postobj['user_email'] != '':
            MailNotification(key).start()
    resultObj = UploadedFile.objects.get(key=key)
    return redirect(resultObj)

accessed = {}

def getResults(request):
    #TODO: show error on UI not related to sanity-error
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    key = getKeyFromUrl(request)
    obj = UploadedFile.objects.get(key=key)

    #errors: when errors from model UploadedFile set to true -> sanity error, when None -> Error in the pipeline

    showResultCSS = ''
    showProgressBar = None
    refresh = None
    showErrorMessage = None
    errorSeverityPU = None
    errorSummaryPU = None
    errorMessagePU = None
    queuePos = None
    queueLen = None

    if obj.finished_bins == obj.total_bins and obj.total_bins != 0:
        try:
            numAccessed = accessed[key]
            print(numAccessed)
        except:
            numAccessed = 0
        numAccessed += 1
        accessed[key] = numAccessed
        print(accessed[key])
        showResultCSS = 'block'
        showProgressBar = False
        refresh = False
        if(obj.errors == None):
            showErrorMessage = True
            errorSeverityPU = 'error'
            errorSummaryPU = 'unknown error'
            errorMessagePU = 'phendb service not working'
        elif(obj.errors == True):
            showErrorMessage = True
            errorSeverityPU = 'warn'
            errorSummaryPU = 'Invalid input file(s)'
            errorMessagePU = 'Please check the invalid_input_files.log file'
    else:
        numAccessed = 0
        showResultCSS = 'none'
        if(obj.errors == False or obj.errors == True):
            showProgressBar = True
            refresh = True
            showErrorMessage = False
        elif (obj.errors == None):
            refresh = False
            showErrorMessage = True
            showProgressBar = False
            errorSeverityPU = 'error'
            errorSummaryPU = 'unknown error'
            errorMessagePU = 'phendb service not working'

    print(obj.finished_bins)
    print(obj.total_bins)

    #position in queue
    queuePos, queueLen = get_current_position(key)
    if queuePos == None:
        showErrorMessage = True
        errorSeverityPU = 'error'
        errorSummaryPU = 'unknown error'
        errorMessagePU = 'phendb service not working'
        refresh = False
        showProgressBar = False


    context = {'result' : 'download/',
               'showResultCSS' : showResultCSS,
               'showNotification' : True if numAccessed == 1 else False,
               'showProgressBar' : showProgressBar,
               'progress' : (obj.finished_bins * 1.0 / obj.total_bins) * 100 if (obj.total_bins!=0) else 0.001,
               'finished_bins' : str(obj.finished_bins),
               'total_bins' : str(obj.total_bins),
               'refresh' : refresh,
               'showInputFormCSS': 'none',
               'showErrorMessage': showErrorMessage,
               'errorSeverityPU' : errorSeverityPU,
               'errorSummaryPU' : errorSummaryPU,
               'errorMessagePU' : errorMessagePU,
               'queuePos' : queuePos + 1,
               'queueLen' : queueLen}

    pprint(context)

    return HttpResponse(template.render(context, request))

def fileDownload(request):
    key = getKeyFromUrl(request)
    resFile = UploadedFile.objects.get(key = key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="{k}.zip"'.format(k=key)
    return response