import json
from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel

class ResultFigures:
    def initialize(job):
        all_bins = BinInJob.objects.filter(job=job)
        db_data = dict(map(lambda x: PicaResult.objects.filter(bin=x.bin), all_bins))
        resultFigures = ResultFigures()
        resultFigures.__buildTraitPredictionTable(db_data)

    def __buildTraitPredictionTable(self, db_data):
        #TODO: change
        json_data = json.dumps(db_data)
        return json_data

