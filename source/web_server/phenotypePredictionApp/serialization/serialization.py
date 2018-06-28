from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel

class PicaResultForUI:

    def __init__(self, job):
        self.job = job
        self.all_bins = BinInJob.objects.filter(job=job)

    def test(self):
        return self.__table_calc_trait_prediction()

    def __table_calc_trait_prediction(self):
        obj = TableCalcTraitPrediction(self)
        return obj.calc()


    def __table_calc_prediction_details(self):
        pass

    def __table_calc_trait_counts(self):
        pass


class TableCalcTraitPrediction:

    titles = {'bin' : 'Bin'}

    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI


    def calc(self):
        arr = []
        for bin_obj in self.picaResultForUI.all_bins:
            single_pica_result = PicaResult.objects.filter(bin=bin_obj.bin)
            bin_name = BinInJob.objects.get(bin=bin_obj.bin, job=self.picaResultForUI.job)
            arr.append(self.__parse_bin(single_pica_result, bin_name))
        return arr


    def __parse_bin(self, single_pica_result, bin_name):
        arr = []
        for item in single_pica_result:
            arr.append("" if bin_name.bin_alias is None else bin_name.bin_alias)
            arr.append(item.model.model_name)
            arr.append("+" if item.verdict else "-")
            arr.append(item.pica_pval)
            arr.append(item.accuracy)
        return arr
