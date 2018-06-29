from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel

class PicaResultForUI:

    def __init__(self, job):
        self.job = job
        self.all_bins_in_job = BinInJob.objects.filter(job=job)
        self.__calc_prediction_details()


    def __calc_prediction_details(self):
        self.prediction_details = _PredictionDetails(self)


    def __calc_prediction(self):
        self.prediction = _Prediction(self)

    def __calc_trait_counts(self):
        pass


class _PredictionDetails:

    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title" : "Bin"}, {"title" : "Model"}, {"title" : "Prediction"}, {"title" : "Prediction-Confidence"}, {"title" : "Balanced-Accuracy"}]

    def get_values(self):
        return self.values

    def get_titles(self):
        return _PredictionDetails.TITLES

    def __calc(self):
        arr = []
        for bin_in_job in self.picaResultForUI.all_bins_in_job:
            single_pica_result = PicaResult.objects.filter(bin=bin_in_job.bin)
            bin_name = bin_in_job.bin_alias
            arr = arr + self.__parse_bin(single_pica_result, bin_name)
        self.values = arr


    def __parse_bin(self, single_pica_result, bin_name):
        arr = []
        for item in single_pica_result:
            single_row = []
            single_row.append("" if bin_name is None else bin_name)
            single_row.append(item.model.model_name)
            single_row.append("+" if item.verdict else "-")
            single_row.append(item.pica_pval)
            single_row.append(item.accuracy)
            arr.append(single_row)
        return arr

class _Prediction:

    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    def __calc(self):
        for bin_in_job in self.picaResultForUI.all_bins_in_job:
            bin_name = bin_in_job.bin_alias
            bin = bin_in_job.bin
            pica_results = PicaResult.objects.filter(bin=bin)
            pica_models = PicaModel.objects.all()
            self.values = [bin_name]
            self.titles = [""]
            for pica_model in pica_models:
                pica_result = PicaResult.objects.filter(bin=bin, model=pica_model)
                if len(pica_result) == 0:
                    continue #model not used in this prediction (e.g. old model)
                self.values.append("+" if pica_result.verdict else "-")
                self.titles.append(pica_result.model.model_name)

    def get_values(self):
        return self.values

    def get_titles(self):
        return self.titles