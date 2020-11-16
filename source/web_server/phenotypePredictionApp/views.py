from django.http import HttpResponse
from django.template import loader
from .forms import FileForm
from .variables import PHENDB_BASEDIR, PHENDB_DATA_DIR, PHENDB_QUEUE, PHENDB_DEBUG
from django.shortcuts import redirect
from businessLogic.mailNotification import MailNotification
import uuid
from businessLogic.startProcess import StartProcessThread
from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel
from ipware.ip import get_real_ip
from redis import Redis
from rq import Queue, get_current_job
from .serialization.serialization import PicaResultForUI
import struct
from phenotypePredictionApp.global_variables import DEFAULT_VALUES
from django.http import JsonResponse
from django.core.files import File
import os



#------------------functions---------------------------------------------
#useful functions, NO views

def getKeyFromUrl(request):
    reqPath = request.path
    if(reqPath[-1] == '/'):
        reqPath = reqPath[:-1]
    urlparts = reqPath.rsplit('/')
    for part in urlparts:
        if len(part) == 36 or part == "PHENDB_PRECALC":
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

def get_sample_fileobj():
    example_file_path = DEFAULT_VALUES["example_file_path"]
    example_file_name = os.path.basename(example_file_path)
    example_file = open(example_file_path, "rb")
    django_example_file = File(example_file)
    return django_example_file, example_file_name


#-------------------Views-------------------------------------------------

def index(request):
    template = loader.get_template('phenotypePredictionApp/index.xhtml')
    context = {'result' : 'No result yet',
               'showResultCSS' : 'none',
               'showNotification' : False,
               'showInputFormCSS': 'block',
               'showProgressBar' : False,
               'refresh' : False,
               'queueLen' : getQueueLength(),
               'requested_balac_default' : DEFAULT_VALUES['balanced_accuracy_cutoff'],
               'requested_conf_default' : DEFAULT_VALUES['prediction_confidence_cutoff'],
               'example_data_alert': DEFAULT_VALUES['example_file_alert'],
               }
    return HttpResponse(template.render(context, request))

def sendinput(request):
    postobj = request.POST.copy()
    fileobj = request.FILES.copy()
    key = str(uuid.uuid4())
    postobj['key'] = key

    if postobj.get('example_data_button', None) is not None:
        example_file, example_name = get_sample_fileobj()
        postobj['filename'] = example_name
        postobj['fileInput'] = example_file
    else:
        # works only if just one file is uploaded
        for filename, file in request.FILES.items():
            name = request.FILES[filename].name
        postobj['filename'] = name  # this will return error if example data is used (if/else logic needed)
        postobj['fileInput'] = fileobj['fileInput']
    postobj['user_ip'] = get_real_ip(request)
    postobj['errors'] = False

    if postobj.get('requested_balac', None) is None:
        postobj['requested_balac'] = DEFAULT_VALUES['balanced_accuracy_cutoff']
        postobj['requested_conf'] = DEFAULT_VALUES['prediction_confidence_cutoff']

    form = FileForm(postobj, fileobj)

    if(form.is_valid()):
        modelInstance = form.save(commit=False)
        modelInstance.save()
        StartProcessThread(key).start()
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
        job = Job.objects.get(key=key)
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

    pica_result_for_ui = None

    if job.finished_bins == job.total_bins and job.total_bins != 0:
        try:
            numAccessed = accessed[key]
        except:
            numAccessed = 0
        numAccessed += 1
        accessed[key] = numAccessed
        showResultCSS = 'block'
        showProgressBar = False
        refresh = False
        if(job.error_type == "UNKNOWN"):
            showErrorMessage = True
            errorSeverityPU = 'error'
            errorSummaryPU = 'Unknown Error'
            errorMessagePU = 'An unknown internal error has occurred.'
        elif(job.error_type == "INPUT"):
            showErrorMessage = True
            errorSeverityPU = 'warn'
            errorSummaryPU = 'Invalid Input File(s)'
            errorMessagePU = 'Please check the invalid_input_files.log file.'
        #results to display in UI
        all_bins = BinInJob.objects.filter(job=job)
        #TODO: resultsList initialzation or other method to get all results for UI
        pica_result_for_ui = PicaResultForUI(job=job)
    else:
        numAccessed = 0
        showResultCSS = 'none'
        if(job.error_type in ("INPUT", "")):
            showProgressBar = True
            refresh = True
            showErrorMessage = False
        elif (job.error_type == "UNKNOWN"):
            refresh = False
            showErrorMessage = True
            showProgressBar = False
            errorSeverityPU = 'error'
            errorSummaryPU = 'Unknown Error'
            errorMessagePU = 'An unknown internal error has occurred.'
        elif (job.error_type == "ALL_DROPPED"):
            refresh = False
            showErrorMessage = True
            showProgressBar = False
            errorSeverityPU = 'error'
            errorSummaryPU = 'No Valid Input Files'
            errorMessagePU = 'None of the uploaded files were valid input files for PhenDB.'


    #position in queue
    queuePos, queueLen = get_current_position(key)

    # write queue length to binary file
    with open(os.path.join(PHENDB_DATA_DIR, "logs/queuelen"), "wb") as bytefile:
        for times in range(queueLen):
            bytefile.write(struct.pack('x'))

    if queuePos == None:
        showErrorMessage = True
        errorSeverityPU = 'error'
        errorSummaryPU = 'Unknown Error'
        errorMessagePU = 'An unknown internal error has occurred.'
        refresh = False
        showProgressBar = False

    context = {'result' : 'download/',
               'showResultCSS' : showResultCSS,
               'showNotification' : True if numAccessed == 1 else False,
               'showProgressBar' : showProgressBar,
               'progress' : (job.finished_bins * 1.0 / job.total_bins) * 100 if job.total_bins != 0 else 0, # necessary to avoid DivBy0 Exception
               'finished_bins' : str(job.finished_bins),
               'total_bins' : str(job.total_bins),
               'refresh' : refresh,
               'showInputFormCSS': 'none',
               'showErrorMessage': showErrorMessage,
               'errorSeverityPU' : errorSeverityPU,
               'errorSummaryPU' : errorSummaryPU,
               'errorMessagePU' : errorMessagePU,
               'queuePos' : queuePos + 1,
               'queueLen' : queueLen,
               'prediction_details_values' : pica_result_for_ui.prediction_details.get_values() if pica_result_for_ui is not None else "",
               'prediction_details_titles' : pica_result_for_ui.prediction_details.get_titles() if pica_result_for_ui is not None else "",
               'prediction_values' : pica_result_for_ui.prediction.get_values() if pica_result_for_ui is not None else "",
               'prediction_titles' : pica_result_for_ui.prediction.get_titles() if pica_result_for_ui is not None else "",
               'trait_counts_values' : pica_result_for_ui.trait_counts.get_values() if pica_result_for_ui is not None else "",
               'trait_counts_titles': pica_result_for_ui.trait_counts.get_titles() if pica_result_for_ui is not None else "",
               'bin_summary_values' : pica_result_for_ui.bin_summary.get_values() if pica_result_for_ui is not None else "",
               'bin_summary_titles': pica_result_for_ui.bin_summary.get_titles() if pica_result_for_ui is not None else "",
               'bin_alias_list' : pica_result_for_ui.bin_alias_list if pica_result_for_ui is not None else "",
               'model_list' : pica_result_for_ui.prediction.get_raw_title_list() if pica_result_for_ui is not None else "",
               }

    return HttpResponse(template.render(context, request))

def updateResultsAjax(request):
    disable_cutoffs = False if request.GET.get('disable_cutoffs') is None else True
    requested_balac = request.GET.get('requested_balac')
    requested_conf = request.GET.get('requested_conf')

    key = getKeyFromUrl(request)
    job = Job.objects.get(key=key)
    pica_result_for_ui = PicaResultForUI(job=job, requested_conf=requested_conf, requested_balac=requested_balac, disable_cutoffs=disable_cutoffs)

    data = {'prediction_details_values': pica_result_for_ui.prediction_details.get_values() if pica_result_for_ui is not None else "",
            'prediction_details_titles': pica_result_for_ui.prediction_details.get_titles() if pica_result_for_ui is not None else "",
            'prediction_values': pica_result_for_ui.prediction.get_values() if pica_result_for_ui is not None else "",
            'prediction_titles': pica_result_for_ui.prediction.get_titles() if pica_result_for_ui is not None else "",
            'trait_counts_values': pica_result_for_ui.trait_counts.get_values() if pica_result_for_ui is not None else "",
            'trait_counts_titles': pica_result_for_ui.trait_counts.get_titles() if pica_result_for_ui is not None else "",
            'bin_summary_values': pica_result_for_ui.bin_summary.get_values() if pica_result_for_ui is not None else "",
            'bin_summary_titles': pica_result_for_ui.bin_summary.get_titles() if pica_result_for_ui is not None else "",
            'bin_alias_list': pica_result_for_ui.bin_alias_list if pica_result_for_ui is not None else "",
            'model_list': pica_result_for_ui.prediction.get_raw_title_list() if pica_result_for_ui is not None else "",
    }
    return JsonResponse(data)

def fileDownload(request):
    key = getKeyFromUrl(request)
    resFile = Job.objects.get(key=key)
    response = HttpResponse(resFile.fileOutput, content_type='application/tar+gzip')
    response['Content-Disposition'] = 'attachment; filename="phendb_{k}.zip"'.format(k=key)
    return response


# ---------------VIEWS-ERRORS--------------------------------------

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
