#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 6/13/18.
#
import django
import datetime
import os
import sys  # remove this in process

sys.path.append("/apps/phenDB/source/web_server")  # remove this in process
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import Job, Bin, BinInJob, PicaModel, PicaResult, Taxon

BALAC_CUTOFF = 0.75  # "${params.accuracy_cutoff}"
PICA_CONF_CUTOFF = 0.75  # "${params.pica_conf_cutoff}"
TRAIT_DEPENDENCY_FILE = "trait_dependencies.tsv"
SHOW_ALL_RESULTS = "false" == "true"  # "${params.show_everything}" == "true"
ROUND_TO = 2
BIN_MDSUMS = ["0ccb1344d608c883164aa949ad1a06bb",
              "a8f3265e27ea5f4af49f03802f12b1ac",
              "5f341600ffec0b8c689cba69d412aa82",
              "d24a691dc232586ce893852c4e3f7ee7",
              "20dc8a09dc90c0db18d67ff2284b035c"]  # sorted("${mdsum_string}".split("\t"), reverse=True)
JOB_ID = "b275cd2e-8aaa-46fa-97f0-d416a3c5e14b" # "${jobname}"

now = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
INDIVIDUAL_RESULTS_HEADER = "Model_Name\tPrediction\tPrediction_Confidence\tBalanced_Accuracy\nModel_Description\n"
BIN_SUMMARY_HEADER = "Bin_Name\tCompleteness\tContamination\tStrain_Heterogeneity\tTaxon_ID\tTaxon_Name\tTaxon_Rank\n"
KRONA_FILE_HEADER = "#bin name\t#taxon_id\n"


def filter_by_cutoff(rd, rl, balac, pica, show_all=False):
    if show_all:
        return rd

    for res in rl:
        md5sum = res.bin.md5sum
        model = res.model.model_name
        if (res.accuracy <= balac or res.pica_pval <= pica) and model != "ARCHAEA":
            if rd[md5sum][model]["PICA_probability"] is not None:
                rd[md5sum][model]["Prediction"] = "n.d."
    return rd


def filter_by_hierarchy(rd, bl, ml, schema, show_all=False):
    # fix taxonomy for bins predicted to be archaea
    for bin in bl:
        if rd[bin.md5sum]["ARCHAEA"]["Prediction"] == "+":
            bin.tax_id = "2157"
            bin.save()
    if show_all:
        return rd, bl

    # read in dependencies
    dep_dic = {model_name: {} for model_name in ml}
    with open(schema, "r") as depfile:
        for line in depfile:
            model_deps = line.strip().split("\t")
            mname = model_deps[0]
            deps = {dep.split(":")[0]: dep.split(":")[1] for dep in model_deps[1:]}
            if mname.startswith("!"):
                mname = mname[1:]
                for name in ml:
                    if name != mname:
                        dep_dic[name] = {**dep_dic[name], **deps}
            else:
                dep_dic[mname] = {**dep_dic[mname], **deps}

    pprint(dep_dic)

    # lookup dependencies for each bin for each model and set to ND if not all satisfied
    for bin in bl:
        for model_name in ml:
            current_cstr = dep_dic[model_name]
            for constraint_model, c_s in current_cstr.items():
                try:
                    if rd[bin.md5sum][constraint_model]["Prediction"] not in ("n.d.", "n.c.", c_s):
                        rd[bin.md5sum][model_name]["Prediction"] = "n.c."
                        # also update all results of this bin with this model (also old ones)
                        # n.c. masking in db, for web UI
                        to_be_masked = PicaResult.objects.filter(bin=bin, model__model_name=model_name)
                        to_be_masked.update(nc_masked=True)
                except KeyError:
                    print("Invalid constraint definition: unknown model.")
    return rd, bl  # remove in the end


# model-bin-related information
parentjob = Job.objects.get(key=JOB_ID)
job_bins = list(Bin.objects.filter(md5sum__in=BIN_MDSUMS))
job_bins_aliases = [BinInJob.objects.get(bin=x, job=parentjob).bin_alias for x in job_bins]
jam = {x.md5sum: y for x, y in zip(job_bins, job_bins_aliases)}
all_job_results = PicaResult.objects.filter(bin__in=job_bins)
first_bin_results = PicaResult.objects.filter(bin=job_bins[0])

# set ARCHAEA to confidence 1 ==> no Platt scaling in current model
archaea_results = all_job_results.filter(model__model_name="ARCHAEA")
archaea_results.update(pica_pval="1")

all_model_names = sorted(list(set([x.model.model_name for x in first_bin_results])))
all_newest_models = [PicaModel.objects.filter(model_name=x).latest("model_train_date") for x in all_model_names]
all_model_desc = [x.model_desc for x in all_newest_models]
job_results_dict = {x.md5sum: {y: None for y in all_model_names} for x in job_bins}
model_results_count = {x: {"+": 0, "-": 0, "n.d.": 0, "n.c.": 0} for x in all_model_names}

for pica_result in all_job_results:
    pica_verdict_to_string = {1: "+", 0: "-"}
    string_verdict = pica_verdict_to_string[int(pica_result.verdict)]
    this_result_dict = {
        "Model_name"       : pica_result.model.model_name,
        "Model_description": pica_result.model.model_desc,
        "Prediction"       : string_verdict,
        "PICA_probability" : pica_result.pica_pval,
        "Balanced_Accuracy": pica_result.accuracy
    }
    job_results_dict[pica_result.bin.md5sum][pica_result.model.model_name] = this_result_dict

# filter results as desired
cutoff_filtered_job_results = filter_by_cutoff(job_results_dict,
                                               all_job_results,
                                               balac=BALAC_CUTOFF,
                                               pica=PICA_CONF_CUTOFF,
                                               show_all=SHOW_ALL_RESULTS)

from pprint import pprint

hf_job_results, bins_tax_fixed = filter_by_hierarchy(cutoff_filtered_job_results,
                                                     job_bins,
                                                     all_model_names,
                                                     schema=TRAIT_DEPENDENCY_FILE,
                                                     show_all=SHOW_ALL_RESULTS)


# write summaries and individual files
with open("trait_summary_matrix.csv", "w") as summary_matrix:
    summary_matrix.write("# PhenDB\n"
                         "# Time of run: {da}\n"
                         "# Accuracy cut-off: {co}\n"
                         "# Display all results? {unl}\n".format(da=now,
                                                                 co=BALAC_CUTOFF,
                                                                 unl=SHOW_ALL_RESULTS))
    summary_matrix.write("# Summary of Trait Prediction Results:\n\n")
    summary_matrix.write("\t" + "\t".join(all_model_desc) + "\n")
    summary_matrix.write("Bin Name\t" + "\t".join(all_model_names) + "\n")

    # write individual result files and append to list for each bin
    for bin in bins_tax_fixed:
        all_bin_predictions = []
        with open("{bn}.traits.csv".format(bn=jam[bin.md5sum]), "w") as ind_t:
            ind_t.write(INDIVIDUAL_RESULTS_HEADER)
            for model in all_model_names:
                result_for_write = hf_job_results[bin.md5sum][model]
                if type(result_for_write["PICA_probability"]) == float:
                    result_for_write["PICA_probability"] = round(result_for_write["PICA_probability"], ROUND_TO)
                if type(result_for_write["Balanced_Accuracy"]) == float:
                    result_for_write["Balanced_Accuracy"] = round(result_for_write["Balanced_Accuracy"], ROUND_TO)
                all_bin_predictions.append(result_for_write["Prediction"])
                model_results_count[model][result_for_write["Prediction"]] += 1
                ind_t.write("{mn}\t{pred}\t{pica_conf}\t{balac}\t{desc}\n".format(mn=result_for_write["Model_name"],
                                                                                  pred=result_for_write["Prediction"],
                                                                                  pica_conf=result_for_write[
                                                                                      "PICA_probability"],
                                                                                  balac=result_for_write[
                                                                                      "Balanced_Accuracy"],
                                                                                  desc=result_for_write[
                                                                                      "Model_description"]
                                                                                  ))
        summary_matrix.write("\t".join([jam[bin.md5sum], *all_bin_predictions]) + "\n")

with open("trait_counts.csv", "w") as count_table:
    count_table.write("Model Name\t+\t-\tn.d.\tn.c.\n")
    for model in all_model_names:
        y = str(model_results_count[model]["+"])
        n = str(model_results_count[model]["-"])
        nd = str(model_results_count[model]["n.d."])
        nc = str(model_results_count[model]["n.c."])
        count_table.write("\t".join([model, y, n, nd, nc]) + "\n")

# summarize bin-related information
with open("bin_summary.csv", "w") as taxfile:
    taxfile.write(BIN_SUMMARY_HEADER)
    for bin in bins_tax_fixed:
        try:
            taxon_entry = Taxon.objects.get(tax_id=bin.tax_id)
            tid = taxon_entry.tax_id
            tname = taxon_entry.taxon_name
            trank = taxon_entry.taxon_rank
        except ObjectDoesNotExist:
            tid = "1"
            tname = "root"
            trank = "no rank"
        taxfile.write("{bn}\t{com}\t{con}\t{sh}\t{ti}\t{tn}\t{tr}\n".format(bn=jam[bin.md5sum],
                                                                            com=bin.comple,
                                                                            con=bin.conta,
                                                                            sh=bin.strainhet,
                                                                            ti=tid,
                                                                            tn=tname,
                                                                            tr=trank))
with open("taxonomy_krona.tsv", "w") as tax_krona:
    tax_krona.write(KRONA_FILE_HEADER)
    for bin in bins_tax_fixed:
        tax_krona.write("{bn}\t{ti}\n".format(bn=jam[bin.md5sum],
                                              ti=bin.tax_id))