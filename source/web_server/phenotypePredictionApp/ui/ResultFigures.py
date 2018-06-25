import json
from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel
from django.core.serializers import serialize

class ResultFigures:
    def initialize(job):
        all_bins = BinInJob.objects.filter(job=job)
        db_data = list(map(lambda x: PicaResult.objects.filter(bin=x.bin), all_bins))
        resultFigures = ResultFigures()
        #resultFigures.__buildTraitPredictionTable(db_data)
        return resultFigures

    def __buildTraitPredictionTable(self, db_data):
        #TODO: change
        self.traitPredictionTable = serialize('json', db_data, use_natural_foreign_keys=True)

    def getTraitPredictionTable(self):
        return self.traitPredictionTable
