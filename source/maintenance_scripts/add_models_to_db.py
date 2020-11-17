#!/usr/bin/env python3
from typing import Union, Optional, Tuple
from pathlib import Path
import json

import click
import django
from django.utils import timezone


django.setup()
from phenotypePredictionApp.models import (
    Enog, PicaModel, EnogRank, PicaModelAccuracy, PicaModelTrainingData
)

DB_ENOGS = {v.enog_name: v for k, v in Enog.objects.in_bulk().items()}


def load_model(d: Union[str, Path]) -> Optional[Tuple]:
    d = Path(str(d))

    with open(d/f'{d.name}.phenotype') as fin_p, open(d/f'{d.name}.taxids') as fin_t:
        accession_type, trait_name = fin_p.readline().strip().split('\t')
        train_data = [x.strip().split() for x in fin_p.readlines()]
        taxid_dict = {
            x: y
            for x, y in [x.strip().split('\t') for x in fin_t.readlines()]
        }
        for k, v in list(taxid_dict.items()):
            gcf = k.replace('GCF_', 'GCA_')
            taxid_dict[gcf] = v

    try:
        with open(d/f'{d.name}.description') as fin:
            desc = fin.read().replace('\n', ' ')
    except FileNotFoundError:
        desc = ''

    try:
        with open(d/f'{d.name}.type') as fin:
            model_type = fin.read().replace('\n', '')
    except FileNotFoundError:
        model_type = 'phenotrex-TrexSVM'

    newmodel = PicaModel(
        model_name=trait_name,
        model_desc=desc,
        alg_type=model_type,
        feature_type='eggNOG5-tax-2',
        model_train_date=timezone.now()
    )
    try:
        newmodel.save()
    except Exception as e:
        print(f'Error during model saving: {e}')
        raise e


    with open(d/f'{d.name}.rank') as fin:
        fin.readline()
        wts = [x.strip().split() for x in fin.readlines()]
    enog_rank_list = [
        EnogRank(
            model=newmodel,
            enog=DB_ENOGS[x],
            internal_rank=i + 1,
            score=float(y),
            pred_class=(False if float(y) <= -1 else True))
        for i, (_, x, y) in enumerate(wts)
    ]
    try:
        EnogRank.objects.bulk_create(enog_rank_list)
    except Exception as e:
        print(f'Error during EnogRank data saving: {e}')
        newmodel.delete()
        raise e

    with open(d/f'{d.name}.accuracy.json') as fin:
        accs = json.load(fin)

    acc_list = [
        PicaModelAccuracy(
            model=newmodel,
            comple=accs[i]['completeness'],
            conta=accs[i]['contamination'],
            mean_balanced_accuracy=accs[i]['mean_balanced_accuracy'],
            mean_fp_rate=None,
            mean_fn_rate=None
        )
        for i in range(21 * 21)
    ]
    try:
        PicaModelAccuracy.objects.bulk_create(acc_list)
    except Exception as e:
        print(f'Error during PicaModelAccuracy saving: {e}')
        newmodel.delete()
        raise e

    train_list = [
        PicaModelTrainingData(
            model=newmodel,
            tax_id=taxid_dict[x],
            assembly_id=x,
            verdict=(True if y == 'YES' else False)
        )
        for x, y in train_data
    ]
    try:
        PicaModelTrainingData.objects.bulk_create(train_list)
    except Exception as e:
        print(f'Error during PicaModelTrainingData saving: {e}')
        newmodel.delete()

    return newmodel, train_list, enog_rank_list, acc_list


@click.command()
@click.argument('model_dir', type=click.Path(file_okay=False, exists=True))
@click.option('--drop', is_flag=True)
def main(model_dir: Union[str, Path], drop):
    if drop:
        print('Deleting all existing models and assorted data from DB .')
        PicaModel.objects.all().delete()

    results = []
    for model in Path(str(model_dir)).iterdir():
        if model.is_file():
            continue
        try:
            print(f'Attempting to upload model {model}...')
            load_model(model)
        except Exception as e:
            print(f'Error during Model data loading: {e}')
            results.append((model, False))
            continue
        results.append((model, True))
    for m, r in results:
        if not r:
            print(f'Upload failed for model {m}.')


if __name__ == '__main__':
    main()