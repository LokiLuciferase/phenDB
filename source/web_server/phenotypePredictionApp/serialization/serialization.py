from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel

class PicaResultForUI:

    def __init__(self, job):
        self.job = job
        self.all_bins = BinInJob.objects.filter(job=job)

    def test(self):
        self.__table_calc_trait_prediction()

    def __table_calc_trait_prediction(self):
        for bin_obj in self.all_bins:
            single_pica_result = PicaResult.objects.filter(bin=bin_obj.bin)
            bin_name = BinInJob.objects.filter(bin=bin_obj.bin, job=self.job)
            test = zip(single_pica_result, bin_name)
            print(test)


    def __table_calc_prediction_details(self):
        pass

    def __table_calc_trait_counts(self):
        pass

    '''def __calc

def picaResultForUI(job):
    resultList = []
    all_bins = BinInJob.objects.filter(job=job)
    for bin_obj in all_bins:
        resultsList += PicaResult.objects.filter(bin=bin_obj.bin).select_related('bin') '''
