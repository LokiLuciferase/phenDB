from phenotypePredictionApp.models import Job, PicaResult, BinInJob, PicaModel, Taxon

class PicaResultForUI:

    def __init__(self, job, requested_conf=None, requested_balac=None, disable_cutoffs=None):
        self.job = job
        self.job_date = self.job.job_date
        self.requested_balac = float(job.requested_balac) if not requested_balac else float(requested_balac)
        self.requested_conf = float(job.requested_conf) if not requested_conf else float(requested_conf)
        self.disable_cutoffs = bool(job.disable_cutoffs) if not disable_cutoffs else bool(disable_cutoffs)
        self.all_bins_in_job = BinInJob.objects.filter(job=job).select_related('bin')
        self.all_bins = [x.bin for x in self.all_bins_in_job]
        self.bin_alias_list = [x.bin_alias for x in self.all_bins_in_job]
        print("fetched all bins")
        all_models_for_currentjob = [x.model for x in PicaResult.objects.filter(bin=self.all_bins[0]).select_related('model')]
        all_model_names_for_cj = set()
        self.newest_models_for_currentjob = []
        # sort by model train date and add newest possible results for models older than job
        for pm in sorted(all_models_for_currentjob, key=lambda x: x.model_train_date, reverse=True):
            pmn = pm.model_name
            pmd = pm.model_train_date
            if pmd > self.job_date:
                continue
            if pmn not in all_model_names_for_cj:
                all_model_names_for_cj.add(pmn)
                self.newest_models_for_currentjob.append(pm)
        self.all_results_for_currentjob = PicaResult.objects.filter(bin__in=self.all_bins, model__in=self.newest_models_for_currentjob)
        self.resdic = self.__make_results_dict()
        self.__calc_prediction_details()
        self.__calc_prediction()
        self.__calc_trait_counts()
        self.__calc_bin_summary()

    def __make_results_dict(self):
        resdic = {}
        bin_id_to_bij_id = {x.bin_id: x.bin_alias for x in self.all_bins_in_job}
        bin_id_to_alias = {x.id: bin_id_to_bij_id[x.id] for x in self.all_bins}
        model_id_to_name = {x.id: x.model_name for x in self.newest_models_for_currentjob}
        for res in list(self.all_results_for_currentjob.values()):
            bin_alias = bin_id_to_alias[res.get("bin_id")]
            model_name = model_id_to_name[res.get("model_id")]
            verdict = res.get("verdict")
            mba = res.get("accuracy")
            pconf = res.get("pica_pval")
            nc_masked = res.get("nc_masked")
            bin_in_resdic = resdic.setdefault(bin_alias, {})
            bin_in_resdic[model_name] = {"verdict": verdict,
                                         "accuracy": mba,
                                         "pica_pval": pconf,
                                         "nc_masked": nc_masked}
        return resdic

    def __calc_prediction_details(self):
        self.prediction_details = _PredictionDetails(self)

    def __calc_prediction(self):
        self.prediction = _Prediction(self)

    def __calc_trait_counts(self):
        self.trait_counts = _TraitCounts(self)

    def __calc_bin_summary(self):
        self.bin_summary = _BinSummary(self)

    def _apply_masks(self, resultdic):
        result_string = "+" if resultdic["verdict"] else "-"
        nc_masked = resultdic["nc_masked"]
        nd_masked = resultdic["pica_pval"] <= self.requested_conf or resultdic["accuracy"] <= self.requested_balac

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

    TITLES = [{"title" : "Bin"}, {"title" : "<a href='http://phendb.org/reports/modeloverview'>link text</a>"}, {"title" : "Prediction"}, {"title" : "Prediction_Confidence"}, {"title" : "Balanced_Accuracy"}]

    def get_values(self):
        return self.values

    def get_titles(self):
        return _PredictionDetails.TITLES

    def __calc(self):
        arr = []
        rd = self.picaResultForUI.resdic
        for bin_name in rd.keys():
            bin_dic = rd[bin_name]
            arr += self.__parse_bin(bin_dic, bin_name)
        self.values = arr

    def __parse_bin(self, bin_dic, bin_name):
        arr = []
        for model_name in sorted(bin_dic.keys()):
            model_dic = bin_dic[model_name]
            single_row = []
            single_row.append(bin_name)
            single_row.append(model_name)
            single_row.append(self.picaResultForUI._apply_masks(model_dic))
            single_row.append(round(model_dic["pica_pval"], 2))
            single_row.append(round(model_dic["accuracy"], 2))
            arr.append(single_row)
        return arr


class _Prediction:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    def __calc(self):
        self.values = []
        rd = self.picaResultForUI.resdic
        for bin_name in rd.keys():
            bin_dic = rd[bin_name]
            values_tmp = [bin_name]
            self.titles = [{"title" : "Bin_name"}]
            self.raw_title_list = []
            for model_name in sorted(bin_dic.keys()):
                model_dic = rd[bin_name][model_name]
                self.titles.append({"title": model_name})
                self.raw_title_list.append(model_name)
                values_tmp.append(self.picaResultForUI._apply_masks(model_dic))
            self.values.append(values_tmp)

    def get_values(self):
        return self.values

    def get_titles(self):
        return self.titles

    def get_raw_title_list(self):
        return self.raw_title_list


class _TraitCounts:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title" : ""}, {"title" : "+"}, {"title" : "-"}, {"title" : "n.d."}, {"title" : "n.c."}]

    def __calc(self):
        self.values = []
        rd = self.picaResultForUI.resdic
        model_names = [x.model_name for x in self.picaResultForUI.newest_models_for_currentjob]
        countdic = {x: {"+":0, "-":0, "n.d.": 0, "n.c.":0} for x in model_names}

        for bin_name in rd.keys():
            bin_dic = rd[bin_name]
            for model_name in bin_dic.keys():
                model_dic = bin_dic[model_name]
                countdic[model_name][self.picaResultForUI._apply_masks(model_dic)] += 1

        for model_name in sorted(model_names):
            plus = countdic[model_name]["+"]
            minus = countdic[model_name]["-"]
            nd = countdic[model_name]["n.d."]
            nc = countdic[model_name]["n.c."]
            self.values.append([model_name, plus, minus, nd, nc])

    def get_values(self):
        return self.values

    def get_titles(self):
        return _TraitCounts.TITLES


class _BinSummary:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title": "Bin"},
              {"title": "Completeness"},
              {"title": "Contamination"},
              {"title": "Strain heterogeneity"},
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
            self.values.append([bin_name, round(comple, 2), round(conta, 2), strainhet, tax_id, taxon_name, taxon_rank])

    def get_values(self):
        return self.values

    def get_titles(self):
        return _BinSummary.TITLES