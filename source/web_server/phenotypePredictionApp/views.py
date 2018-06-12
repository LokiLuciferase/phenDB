from django.shortcuts import render
from django.http import HttpResponse, HttpResponseRedirect
from django.template import loader
from django.urls import reverse
from .forms import FileForm
from .variables import PHENDB_BASEDIR, PHENDB_QUEUE, PHENDB_DEBUG
from django.shortcuts import redirect
from businessLogic.mailNotification import MailNotification
import uuid
from django.core.urlresolvers import resolve
from businessLogic.startProcess import StartProcessThread
from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel
from pprint import pprint
from ipware.ip import get_real_ip
from redis import Redis
from rq import Queue, get_current_job
import struct
import traceback
import os


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
        que = Queue(PHENDB_QUEUE, connection=Redis())
        jobs = list(que.jobs)
        for index, job in enumerate(jobs):
            if job.id == keyname:
                return (index, len(jobs))
        return (-1, len(jobs))
    except:
        return (None, None)

def getQueueLength():
    que = Queue(PHENDB_QUEUE, connection=Redis())
    jobs = list(que.jobs)
    return len(jobs)


#-------------------Views-------------------------------------------------

def index(request):
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResultCSS' : 'none',
               'showNotification' : False,
               'showInputFormCSS': 'block',
               'showProgressBar' : False,
               'refresh' : False,
               'queueLen' : getQueueLength()}
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
        StartProcessThread(key, postobj['requested_balac']).start()
        if postobj['user_email'] != '':
            MailNotification(key).start()
    resultObj = Job.objects.get(key=key)
    return redirect(resultObj)

accessed = {}

def getResults(request):
    #TODO: show error on UI not related to sanity-error
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    templateError = loader.get_template('phenotypePredictionApp/error.xhtml')
    key = getKeyFromUrl(request)
    try:
        obj = Job.objects.get(key=key)
    except Job.DoesNotExist:
        context = {'errorMessage' : 'The requested page does not exist. Please note that all results are deleted after 30 days!'}
        return HttpResponse(templateError.render(context, request))

    #errors: when errors from model Job set to true -> sanity error, when None -> Error in the pipeline

    showResultCSS = ''
    showProgressBar = None
    refresh = None
    showErrorMessage = None
    errorSeverityPU = None
    errorSummaryPU = None
    errorMessagePU = None
    queuePos = None
    queueLen = None
    resultsList = []

    if obj.finished_bins == obj.total_bins and obj.total_bins != 0:
        try:
            numAccessed = accessed[key]
        except:
            numAccessed = 0
        numAccessed += 1
        accessed[key] = numAccessed
        showResultCSS = 'block'
        showProgressBar = False
        refresh = False
        if(obj.error_type == "UNKNOWN"):
            showErrorMessage = True
            errorSeverityPU = 'error'
            errorSummaryPU = 'Unknown Error'
            errorMessagePU = 'An unknown internal error has occurred.'
        elif(obj.error_type == "INPUT"):
            showErrorMessage = True
            errorSeverityPU = 'warn'
            errorSummaryPU = 'Invalid Input File(s)'
            errorMessagePU = 'Please check the invalid_input_files.log file.'
        #results to display in UI
        all_bins = BinInJob.objects.filter(job=obj)
        for bin_obj in all_bins:
            resultsList += PicaResult.objects.filter(bin=bin_obj.bin)
    else:
        numAccessed = 0
        showResultCSS = 'none'
        if(obj.error_type in ("INPUT", "")):
            showProgressBar = True
            refresh = True
            showErrorMessage = False
        elif (obj.error_type == "UNKNOWN"):
            refresh = False
            showErrorMessage = True
            showProgressBar = False
            errorSeverityPU = 'error'
            errorSummaryPU = 'Unknown Error'
            errorMessagePU = 'An unknown internal error has occurred.'
        elif (obj.error_type == "ALL_DROPPED"):
            refresh = False
            showErrorMessage = True
            showProgressBar = False
            errorSeverityPU = 'error'
            errorSummaryPU = 'No Valid Input Files'
            errorMessagePU = 'None of the uploaded files were valid input files for PhenDB.'


    #position in queue
    queuePos, queueLen = get_current_position(key)

    # write queue length to binary file
    with open(os.path.join(PHENDB_BASEDIR, "logs/queuelen"), "wb") as bytefile:
        for times in range(queueLen):
            bytefile.write(struct.pack('x'))

    if queuePos == None:
        showErrorMessage = True
        errorSeverityPU = 'error'
        errorSummaryPU = 'Unknown Error'
        errorMessagePU = 'An unknown internal error has occurred.'
        refresh = False
        showProgressBar = False

    print("len resultsDic " + str(len(resultsList)))
    context = {'result' : 'download/',
               'showResultCSS' : showResultCSS,
               'showNotification' : True if numAccessed == 1 else False,
               'showProgressBar' : showProgressBar,
               'progress' : (obj.finished_bins * 1.0 / obj.total_bins) * 100 if (obj.total_bins!=0 and obj.finished_bins != 0) else 0.001,
               'finished_bins' : str(obj.finished_bins),
               'total_bins' : str(obj.total_bins),
               'refresh' : refresh,
               'showInputFormCSS': 'none',
               'showErrorMessage': showErrorMessage,
               'errorSeverityPU' : errorSeverityPU,
               'errorSummaryPU' : errorSummaryPU,
               'errorMessagePU' : errorMessagePU,
               'queuePos' : queuePos + 1,
               'queueLen' : queueLen,
               'resultsList' : resultsList,
               'all_models' : PicaModel.objects.all}

    #pprint(context)

    return HttpResponse(template.render(context, request))

def fileDownload(request):
    key = getKeyFromUrl(request)
    resFile = Job.objects.get(key=key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="{k}.zip"'.format(k=key)
    return response


#---------------VIEWS-ERRORS--------------------------------------

def permissionDenied(request):
    context = {'errorMessage' : 'You do not have permission to access this site'}
    errorTemplate = loader.get_template('phenotypePredictionApp/error.xhtml')
    return HttpResponse(errorTemplate.render(context,request))

def pageNotFound(request):
    context = {'errorMessage' : 'The requested page does not exist'}
    errorTemplate = loader.get_template('phenotypePredictionApp/error.xhtml')
    return HttpResponse(errorTemplate.render(context,request))

def serverError(request):
    context = {'errorMessage': 'Unknown server error: This should not have happened. Please try again later or contact us!'}
    errorTemplate = loader.get_template('phenotypePredictionApp/error.xhtml')
    return HttpResponse(errorTemplate.render(context, request))
