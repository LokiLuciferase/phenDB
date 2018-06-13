#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 6/13/18.
#
import django
import datetime
import os
import sys  # remove this in process
sys.path.append("/apps/phenDB_devel_LL/source/web_server") # remove this in process
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import Job, Bin, PicaModel, PicaResult


BALAC_CUTOFF = 0.7  # "${params.accuracy_cutoff}"
SHOW_ALL_RESULTS = True if "true" == "true" else False  # add nextflow variable
ROUND_TO = 2
now = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
INDIVIDUAL_RESULTS_HEADER = "Model_name\tPrediction\tBalanced_accuracy\nModel_description\n"
TAXONOMY_HEADER           = "Bin_Name\tTaxon_ID\tTaxon_Name\tTaxon_Rank\n"
KRONA_FILE_HEADER         = "#bin name\t#taxon_id\n"
COMPLECONTA_FILE_HEADER   = "Bin_Name\tCompleteness\tContamination\tStrain_Heterogeneity\n"
BIN_MDSUMS = ["9aedaab39d7a6f59c2e74e8b05c36315",
              "f66b68e21373ea2dbedfa43037178c48",
              "0ccb1344d608c883164aa949ad1a06bb",
              "d24a691dc232586ce893852c4e3f7ee7"]           # sorted("${mdsum_string}".split("\t"), reverse=True)
job_bins = sorted(Bin.objects.filter(md5sum__in=BIN_MDSUMS), key=lambda x: x.bin_name)

# bin-related information
with open("taxonomy.csv", "w") as taxfile:
    taxfile.write(TAXONOMY_HEADER)
    for bin in job_bins:
     taxfile.write("{bn}\t{ti}\t{tn}\t{tr}\n".format(bn=bin.bin_name,
                                                     ti=bin.tax_id,
                                                     tn=bin.taxon_name,
                                                     tr=bin.taxon_rank))

with open("taxonomy_krona.tsv", "w") as tax_krona:
    tax_krona.write(KRONA_FILE_HEADER)
    for bin in job_bins:
     tax_krona.write("{bn}\t{ti}\n".format(bn=bin.bin_name,
                                           ti=bin.tax_id))

with open("completeness_contamination.csv", "w") as cc_out:
    cc_out.write(COMPLECONTA_FILE_HEADER)
    for bin in job_bins:
     cc_out.write("{bn}\t{com}\t{con}\t{sh}\n".format(bn=bin.bin_name,
                                                      com=bin.comple,
                                                      con=bin.conta,
                                                      sh=bin.strainhet))

all_job_results = PicaResult.objects.filter(bin__in=job_bins)
# filter job results here


# model-related information
job_results = all_job_results
all_model_names = sorted(list(set([x.model_name for x in PicaModel.objects.filter()])))
job_results_dict = {x.bin_name: {y: None for y in all_model_names} for x in job_bins}
model_results_count = {x: {"YES": 0, "NO": 0, "ND": 0, "NC": 0} for x in all_model_names}

for pica_result in job_results:
    pica_verdict_to_string = {1: "YES", 0: "NO", 2: "ND", "NC": 3}
    string_verdict = pica_verdict_to_string[int(pica_result.verdict)]
    model_results_count[pica_result.model.model_name][string_verdict] += 1

    this_result_dict = {
        "Model_name": pica_result.model.model_name,
        "Model_description": pica_result.model.model_desc,
        "Prediction": string_verdict,
        "PICA_probability": pica_result.pica_pval,
        "Balanced_Accuracy": pica_result.accuracy
    }
    job_results_dict[pica_result.bin.bin_name][pica_result.model.model_name] = this_result_dict

with open("summary_matrix.csv", "w") as summary_matrix:
    summary_matrix.write("# PhenDB\n"
                         "# Time of run: {da}\n"
                         "# Accuracy cut-off: {co}\n"
                         "# Display unlikely bin-trait combinations? {unl}\n".format(da=now,
                                                                                   co=BALAC_CUTOFF,
                                                                                   unl=SHOW_ALL_RESULTS))
    summary_matrix.write("# Summary of Trait Prediction Results:\n\n")
    summary_matrix.write("Bin Name\t" + "\t".join(all_model_names) + "\n")
    # write individual result files and append to list for each bin
    for bin in job_bins:
        all_bin_predictions = []
        with open("{bn}.traits.csv".format(bn=bin.bin_name), "w") as ind_t:
            ind_t.write(INDIVIDUAL_RESULTS_HEADER)
            for model in all_model_names:
                result_for_write = job_results_dict[bin.bin_name][model]
                if result_for_write is None:
                    all_bin_predictions.append("NC")
                    result_for_write = {
                        "Model_name"       : model,
                        "Model_description": PicaModel.objects.latest(model_name=model).model_desc,
                        "Prediction"       : "NC",
                        "PICA_probability" : "NC",
                        "Balanced_Accuracy": "NC"
                }
                else:
                    all_bin_predictions.append(result_for_write["Prediction"])
                ind_t.write("{mn}\t{pred}\t{balac}\t{desc}\n".format(mn=result_for_write["Model_name"],
                                                                     pred=result_for_write["Prediction"],
                                                                     balac=result_for_write["Balanced_Accuracy"],
                                                                     desc=result_for_write["Model_description"]
                                                                     ))
        summary_matrix.write("\t".join([bin.bin_name, *all_bin_predictions]) + "\n")

with open("trait_counts.csv", "w") as count_table:
    count_table.write("Model Name\tYES\tNO\tND\tNC\n")
    for model in all_model_names:
        y = str(model_results_count[model]["YES"])
        n = str(model_results_count[model]["NO"])
        nd = str(model_results_count[model]["ND"])
        nc = str(model_results_count[model]["NC"])
        count_table.write("\t".join([model, y, n, nd, nc]) + "\n")
