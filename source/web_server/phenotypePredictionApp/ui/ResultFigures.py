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
            resultsDir[bin_obj.bin.bin_name] = PicaResult.objects.filter(bin=bin_obj.bin).values()
        return resultsDir

    def __buildTraitPredictionTable(self, resultsDir):
        #TODO: change
        print(resultsDir)
        #self.traitPredictionTable = serialize('json', resultsDir)
        self.traitPredictionTable = json.dumps(resultsDir)

    def getTraitPredictionTable(self):
        return self.traitPredictionTable
