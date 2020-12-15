from datetime import datetime

import pandas as pd
import click
import django


django.setup()
from phenotypePredictionApp.models import (
    Job,
    BinInJob,
    Bin,
    Enog,
    PicaModel,
    PicaResult,
    PicaResultExplanation,
    PicaModelAccuracy,
    Taxon
)

DB_ENOGS = {v.enog_name: v for k, v in Enog.objects.in_bulk().items()}


def round_5perc(n):
    return round(n * 20) / 20


def check_add_precalc_job(jobname="PHENDB_PRECALC"):
    print("Checking if PHENDB_PRECALC job exists in database...")
    try:
        pc_job = Job.objects.get(key=jobname)
        pc_job.job_date = datetime.now()
        pc_job.save()
        print("PHENDB_PRECALC job date updated.")
    except:
        pc_job = Job(key=jobname)
        pc_job.save()
        print("PHENDB_PRECALC job created.")
    return pc_job


@click.command()
@click.option('--preds')
@click.option('--expl')
@click.option('--ccs')
def import_to_db(preds, expl, ccs):
    preds = pd.read_feather(preds).dropna()
    expl = pd.read_feather(expl).dropna()
    ccs = pd.read_feather(ccs).dropna()
    preds = preds.set_index('Refseq ID')
    expl = expl.set_index('Refseq ID')
    preds['MD5 Sum'] = ccs.set_index('Refseq ID')['MD5 Sum']
    expl['MD5 Sum'] = ccs.set_index('Refseq ID')['MD5 Sum']
    preds.reset_index()
    expl.reset_index()

    singleton_job = check_add_precalc_job()
    print('Creating Bins.')
    all_bins = [
        Bin(
            md5sum=x['MD5 Sum'],
            tax_id=x['Taxon ID'],
            assembly_id=x['Refseq ID'],
            comple=x['Completeness'],
            conta=x['Contamination'],
            strainhet=x['Strain Heterogeneity'],
        )
        for i, x in ccs.iterrows()
    ]
    print('Saving Bins.')
    Bin.objects.bulk_create(all_bins)

    print('Creating BinInJobs.')
    all_bij = [
        BinInJob(
            bin=Bin.objects.get(md5sum=x['MD5 Sum']),
            job=singleton_job,
        )
        for i, x in ccs.iterrows()
    ]
    print('Saving BinsInJobs.')
    BinInJob.objects.bulk_create(all_bij)

    for model_name, df in preds.groupby('Model Name'):
        this_model = PicaModel.objects.filter(model_name=model_name).latest("model_train_date")
        these_preds = []
        print(f'Processing preds for {model_name}...')
        for i, row in df.iterrows():
            bin = Bin.objects.get(md5sum=row['MD5 Sum'])
            comple = round_5perc(bin.comple)
            conta = round_5perc(bin.conta)
            try:
                acc = PicaModelAccuracy.objects.get(
                    model=this_model, comple=comple, conta=conta
                ).mean_balanced_accuracy
                this_pred = PicaResult(
                        bin=bin,
                        model=this_model,
                        verdict=row['Trait present'],
                        pica_pval=row['Confidence'],
                        accuracy=acc,
                )
                these_preds.append(this_pred)
            except:
                print(f'Could not add: bin={bin.assembly_id}, comple={comple}, conta={conta}')
        PicaResult.objects.bulk_create(these_preds)

    for model_name, df in expl.groupby('Model Name'):
        this_model = PicaModel.objects.filter(model_name=model_name).latest("model_train_date")
        these_expl = []
        print(f'Processing expl for {model_name}...')
        for md5sum, ddf in df.groupby('MD5 Sum'):
            this_bin = Bin.objects.get(md5sum=md5sum)
            for j, row in ddf.iterrows():
                this_expl = PicaResultExplanation(
                    bin=this_bin,
                    model=this_model,
                    enog=DB_ENOGS[row['ENOG']],
                    enog_is_present=row['ENOG present'],
                    delta_shap=row['SHAP value']
                )
                these_expl.append(this_expl)
        PicaResultExplanation.objects.bulk_create(these_expl)

if __name__ == '__main__':
    import_to_db()
