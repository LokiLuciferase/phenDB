import json
from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel
from django.core.serializers import serialize

class ResultFigures:
    def initialize(job):
        resultFigures = ResultFigures()
        resultsDir = resultFigures.__createPicaResultDic(job=job)
        resultFigures.__buildTraitPredictionTable(resultsDir)
        return resultFigures

    def __createPicaResultDic(self, job):
        resultsDir = {}
        all_bins = BinInJob.objects.filter(job=job)
        for bin_obj in all_bins:
            resultsDir[bin_obj.bin] = PicaResult.objects.filter(bin=bin_obj.bin)
        return resultsDir

    def __buildTraitPredictionTable(self, resultsDir):
        #TODO: change
        self.traitPredictionTable = serialize('json', resultsDir, use_natural_foreign_keys=True)

    def getTraitPredictionTable(self):
        return self.traitPredictionTable
