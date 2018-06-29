from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel, Taxon

class PicaResultForUI:
    def __init__(self, job, requested_balac=None, requested_conf=None, disable_cutoffs=None):
        self.job = job
        self.requested_balac = float(job.requested_balac) if not requested_balac else float(requested_balac)
        self.requested_conf = float(job.requested_conf) if not requested_conf else float(requested_conf)
        self.disable_cutoffs = bool(job.disable_cutoffs) if not disable_cutoffs else bool(disable_cutoffs)
        self.all_bins_in_job = BinInJob.objects.filter(job=job)
        self.__calc_prediction_details()
        self.__calc_prediction()
        self.__calc_bin_summary()

    def __calc_prediction_details(self):
        self.prediction_details = _PredictionDetails(self)

    def __calc_prediction(self):
        self.prediction = _Prediction(self)

    def __calc_trait_counts(self):
        pass

    def __calc_bin_summary(self):
        self.bin_summary = _BinSummary(self)

    def _apply_masks(self, result):
        result_string = "+" if result.verdict else "-"
        nc_masked = result.nc_masked
        nd_masked = result.pica_pval <= self.requested_conf or result.accuracy <= self.requested_balac

        if self.disable_cutoffs or (not nd_masked and not nc_masked):
            return result_string
        elif nc_masked:
            return "n.c."
        elif nd_masked:
            return "n.d."
        else:
            raise RuntimeError


class _PredictionDetails:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title" : "Bin"}, {"title" : "Model"}, {"title" : "Prediction"}, {"title" : "Pred. Confidence"}, {"title" : "Bal. Accuracy"}]

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
            single_row.append(self.picaResultForUI._apply_masks(item))
            single_row.append(round(item.pica_pval, 2))
            single_row.append(round(item.accuracy, 2))
            arr.append(single_row)
        return arr


class _Prediction:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    def __calc(self):
        self.values = []
        for bin_in_job in self.picaResultForUI.all_bins_in_job:
            bin_name = bin_in_job.bin_alias
            bin = bin_in_job.bin
            pica_models = PicaModel.objects.all()
            values_tmp = [bin_name]
            self.titles = [""]
            for pica_model in pica_models:
                pica_result = PicaResult.objects.filter(bin=bin, model=pica_model)
                if len(pica_result) == 0:
                    continue # model not used in this prediction (e.g. old model)
                self.titles.append({"title": pica_result[0].model.model_name})
                values_tmp.append(self.picaResultForUI._apply_masks(pica_result[0]))

            self.values.append(values_tmp)

    def get_values(self):
        return self.values

    def get_titles(self):
        return self.titles


class _TraitCounts:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title" : ""}, {"title" : "+"}, {"title" : "-"}, {"title" : "n.d."}, {"title" : "n.c."}]

    def __calc(self):
        self.values = []
        pica_models = PicaModel.objects.all()
        for pica_model in pica_models:
            tmp_arr = []
            tmp_arr.append(pica_model.model_name)
            pica_results = PicaResult.objects.filter(model=pica_model)
            if(len(pica_results) == 0):
                continue #model not used in this prediction (e.g. old model)


class _BinSummary:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title": "Bin"},
              {"title": "Completeness"},
              {"title": "Contamination"},
              {"title": "Strain\nhet."},
              {"title": "Taxon ID"},
              {"title": "Taxon name"},
              {"title": "Taxon rank"}]

    def __calc(self):
        self.values = []
        for bin_in_job in self.picaResultForUI.all_bins_in_job:
            bin = bin_in_job.bin
            tax_id = bin.tax_id
            bin_name = bin_in_job.bin_alias
            comple = bin.comple
            conta = bin.conta
            strainhet = bin.strainhet
            taxon = Taxon.objects.get(tax_id=tax_id)
            taxon_name = taxon.taxon_name
            taxon_rank = taxon.taxon_rank
            self.values.append([bin_name, comple, conta, strainhet, tax_id, taxon_name, taxon_rank])

    def get_values(self):
        return self.values

    def get_titles(self):
        return _BinSummary.TITLES