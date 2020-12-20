from phenotypePredictionApp.models import (
    Job, PicaResult, BinInJob, PicaModel, PicaResultExplanation, Taxon, Enog
)


class PicaResultForUI:
    def __init__(self, job, requested_conf=None, requested_balac=None, disable_cutoffs=None):
        self.job = job
        self.job_date = self.job.job_date
        self.requested_balac = (
            float(job.requested_balac) if requested_balac is None else float(requested_balac)
        )
        self.requested_conf = (
            float(job.requested_conf) if requested_conf is None else float(requested_conf)
        )
        self.disable_cutoffs = (
            bool(job.disable_cutoffs) if disable_cutoffs is None else bool(disable_cutoffs)
        )
        self.get_explanations = bool(job.get_explanations)
        self.all_bins_in_job = BinInJob.objects.filter(job=job).select_related("bin")
        self.all_enogs = Enog.objects.all()
        self.all_bins = [x.bin for x in self.all_bins_in_job]
        self.bin_alias_list = [x.bin_alias for x in self.all_bins_in_job]
        all_models_for_currentjob = [
            x.model for x in PicaResult.objects.filter(bin=self.all_bins[0]).select_related("model")
        ]
        if self.get_explanations:
            self.all_expl_in_job = PicaResultExplanation.objects.filter(
                bin__in=self.all_bins, model__in=all_models_for_currentjob
            )
        else:
            self.all_expl_in_job = {}
        all_model_names_for_currentjob = set()
        self.newest_models_for_currentjob = []
        # sort by model train date and add newest possible results for models older than job. Only add unique models
        for picamodel in sorted(
            all_models_for_currentjob, key=lambda x: x.model_train_date, reverse=True
        ):
            picamodel_name = picamodel.model_name
            picamodel_train_date = picamodel.model_train_date
            if picamodel_train_date > self.job_date:
                continue
            if picamodel_name not in all_model_names_for_currentjob:
                all_model_names_for_currentjob.add(picamodel_name)
                self.newest_models_for_currentjob.append(picamodel)
        self.all_results_for_currentjob = PicaResult.objects.filter(
            bin__in=self.all_bins, model__in=self.newest_models_for_currentjob
        )
        self.resdic = self.__make_results_dict()
        self.__calc_prediction_details()
        self.__calc_explanation_details()
        self.__calc_prediction()
        self.__calc_trait_counts()
        self.__calc_bin_summary()

    def __make_results_dict(self):
        """Make a dictionary for each bin represented in the job (via BinInJob)
        and for each PicaModel used in the Job (determined by model_train_date < job_date).
        This reduces database access massively and speeds up loading of web results."""
        resdic = {}
        bin_id_to_bij_id = {x.bin_id: x.bin_alias for x in self.all_bins_in_job}
        bin_id_to_alias = {x.id: bin_id_to_bij_id[x.id] for x in self.all_bins}
        enog_id_to_enog = {x.id: (x.enog_name, x.enog_descr) for x in self.all_enogs}
        model_id_to_name = {x.id: x.model_name for x in self.newest_models_for_currentjob}
        for res in list(self.all_results_for_currentjob.values()):
            model_id = res.get("model_id")
            bin_alias = bin_id_to_alias[res.get("bin_id")]
            model_name = model_id_to_name[model_id]
            verdict = res.get("verdict")
            mean_balac = res.get("accuracy")
            pred_conf = res.get("pica_pval")
            nc_masked = res.get("nc_masked")
            bin_in_resdic = resdic.setdefault(bin_alias, {})

            bin_in_resdic[model_name] = {
                "verdict": verdict,
                "accuracy": mean_balac,
                "pica_pval": pred_conf,
                "nc_masked": nc_masked,
                "model_id": model_id,
                "explanations": []
            }
        for expl in list(self.all_expl_in_job.values()):
            model_id = expl.get('model_id')
            bin_alias = bin_id_to_alias[expl.get("bin_id")]
            model_name = model_id_to_name[model_id]
            enog_id = expl.get('enog_id')
            enog_is_present = expl.get('enog_is_present')
            delta_shap = expl.get('delta_shap')
            enog_name, enog_descr = enog_id_to_enog[enog_id]
            resdic[bin_alias][model_name]['explanations'].append({
                'enog_name': enog_name,
                'enog_descr': enog_descr,
                'enog_is_present': enog_is_present,
                'delta_shap': delta_shap
            })
        return resdic

    def __calc_prediction_details(self):
        self.prediction_details = _PredictionDetails(self)

    def __calc_explanation_details(self):
        self.explanation_details = _ExplanationDetails(self)

    def __calc_prediction(self):
        self.prediction = _Prediction(self)

    def __calc_trait_counts(self):
        self.trait_counts = _TraitCounts(self)

    def __calc_bin_summary(self):
        self.bin_summary = _BinSummary(self)

    def _apply_masks(self, resultdic):
        result_string = "+" if resultdic["verdict"] else "-"
        nc_masked = resultdic["nc_masked"]
        nd_masked = (
            resultdic["pica_pval"] <= self.requested_conf
            or resultdic["accuracy"] <= self.requested_balac
        )

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

    TITLES = [
        {"title": "Bin"},
        {"title": "<a target='_blank' href='http://phendb.org/reports/modeloverview'>Model</a>"},
        {"title": "Prediction"},
        {"title": "Prediction_Confidence"},
        {"title": "Balanced_Accuracy"},
    ]

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
            model_id = model_dic["model_id"]
            single_row = []
            single_row.append(bin_name)
            single_row.append(ModelLink.createModelLink(model_id=model_id, model_name=model_name))
            single_row.append(self.picaResultForUI._apply_masks(model_dic))
            single_row.append(round(model_dic["pica_pval"], 2))
            single_row.append(round(model_dic["accuracy"], 2))
            arr.append(single_row)
        return arr


class _ExplanationDetails:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [
        {"title": "Bin"},
        {"title": "<a target='_blank' href='http://phendb.org/reports/modeloverview'>Model</a>"},
        {"title": "ENOG ID"},
        {"title": "ENOG Present in Bin"},
        {"title": "Resulting SHAP Value"},
    ]

    def get_values(self):
        return self.values

    def get_titles(self):
        return _ExplanationDetails.TITLES

    def __calc(self):
        arr = []
        rd = self.picaResultForUI.resdic
        for bin_name in rd.keys():
            bin_dic = rd[bin_name]
            for model_name in sorted(bin_dic.keys()):
                expls = rd[bin_name][model_name]['explanations']
                for e in expls:
                    arr.append([bin_name, model_name] + list(e.values()))
        self.values = arr


class _Prediction:
    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    def __calc(self):
        self.values = []
        result_dic = self.picaResultForUI.resdic

        for bin_name in result_dic.keys():
            bin_dic = result_dic[bin_name]
            values_tmp = [bin_name]
            self.titles = [{"title": "Bin_name"}]
            self.raw_title_list = []

            for model_name in sorted(bin_dic.keys()):
                model_dic = result_dic[bin_name][model_name]
                model_id = model_dic["model_id"]
                link_with_title = ModelLink.createModelLink(model_id, model_name)
                self.titles.append({"title": link_with_title})
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
    """iterate over all results in result_dict, and increment counters in countdic."""

    def __init__(self, picaResultForUI):
        self.picaResultForUI = picaResultForUI
        self.__calc()

    TITLES = [{"title": ""}, {"title": "+"}, {"title": "-"}, {"title": "n.d."}, {"title": "n.c."}]

    def __calc(self):
        self.values = []
        rd = self.picaResultForUI.resdic
        model_names = [x.model_name for x in self.picaResultForUI.newest_models_for_currentjob]
        countdic = {x: {"+": 0, "-": 0, "n.d.": 0, "n.c.": 0} for x in model_names}

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

    TITLES = [
        {"title": "Bin"},
        {"title": "Completeness"},
        {"title": "Contamination"},
        {"title": "Strain heterogeneity"},
        {"title": "Taxon ID"},
        {"title": "Taxon name"},
        {"title": "Taxon rank"},
    ]

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
            self.values.append(
                [
                    bin_name,
                    round(comple, 2),
                    round(conta, 2),
                    strainhet,
                    tax_id,
                    taxon_name,
                    taxon_rank,
                ]
            )

    def get_values(self):
        return self.values

    def get_titles(self):
        return _BinSummary.TITLES


class ModelLink:
    URL_SKELETON = "http://phendb.org/reports/modeldetails?model_id="

    @staticmethod
    def createModelLink(model_id, model_name):
        link = ModelLink.__getUrlForModel(model_id)
        return "<a target='_blank' href='" + link + "'>" + model_name + "</a>"

    @staticmethod
    def __getUrlForModel(model_id):
        return ModelLink.URL_SKELETON + str(model_id)


class EnogLink:
    URL_SKELETON = "https://phendb.csb.univie.ac.at/reports/enog?enog_name="

    @staticmethod
    def createEnogLink(enog_id):
        link = EnogLink.__getUrlForEnog(enog_id)
        return "<a target='_blank' href='" + link + "'>" + enog_id + "</a>"

    @staticmethod
    def __getUrlForEnog(enog_id):
        return EnogLink.URL_SKELETON + str(enog_id)
#
#
# if __name__ == '__main__':
#
#     from pprint import pprint
#     import django
#     django.setup()
#     from phenotypePredictionApp.models import (
#         Job, PicaResult, BinInJob, PicaModel, PicaResultExplanation, Taxon, Enog
#     )
#     p = PicaResultForUI(job=Job.objects.get(key='ae1fdbff-d331-4240-8b2c-1761c23c760f'))
#     pprint(p.resdic)