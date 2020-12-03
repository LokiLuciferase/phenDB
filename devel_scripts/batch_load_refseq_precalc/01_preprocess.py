#!/usr/bin/env python3
from pathlib import Path
import re

import pandas as pd
import click


CC_FILENAME = 'compleconta_refseq.txt'
TAX_FILENAME = 'taxonomy_mappings.csv'
MD5SUM_FILE = 'md5sums.tsv'
NON_MODEL_FILES = (CC_FILENAME, TAX_FILENAME, MD5SUM_FILE)
PREDS_SUFFIX = '.svm.class.predictions.txt'
EXPL_SUFFIX = '.svm.class.explanations.txt'
REFSEQ_PAT = re.compile(r'GC[AF]_\d+\.\d+')


def load_data(input_directory):
    input_directory = Path(str(input_directory))
    cc_file = input_directory / CC_FILENAME
    model_names = sorted(
        list(
            {
                x.name.split('.')[0] for x in input_directory.iterdir() if x.name not in NON_MODEL_FILES
            }
        )
    )
    preds = []
    expl = []
    for name in model_names:
        preds_df = pd.read_csv(input_directory / f'{name}{PREDS_SUFFIX}', sep='\t', comment='#')
        expl_df = pd.read_csv(input_directory / f'{name}{EXPL_SUFFIX}', sep='\t')
        preds_df['Model Name'] = name
        expl_df['Model Name'] = name
        preds.append(preds_df)
        expl.append(expl_df)
    preds = pd.concat(preds, axis=0)
    expl = pd.concat(expl, axis=0)
    preds.columns = ['Identifier', 'Trait present', 'Confidence', 'Model Name']
    preds['Trait present'] = preds['Trait present'].apply(lambda x: {'YES': True, 'NO': False}[x])

    expl.columns = ['Rank', 'Identifier', 'ENOG', 'ENOG present', 'SHAP value', 'Annotation',
                    'Model Name']
    expl['ENOG present'] = expl['ENOG present'].astype(bool)

    expl.pop('Annotation')
    ccs = pd.read_csv(cc_file, sep='\t', header=None)
    ccs.columns = [
        'Identifier',
        'Completeness',
        'Contamination',
        'Strain Heterogeneity',
        'Found Taxon ID',
        'Found Taxon Name',
        'Found Taxonomic Height'
    ]
    filename_id_mapper = {x: re.findall(REFSEQ_PAT, x)[0] for x in ccs['Identifier']}
    preds['Refseq ID'] = preds['Identifier'].apply(lambda x: filename_id_mapper[x])
    expl['Refseq ID'] = expl['Identifier'].apply(lambda x: filename_id_mapper[x])
    ccs['Refseq ID'] = ccs['Identifier'].apply(lambda x: filename_id_mapper[x])
    tax = pd.read_csv(input_directory / TAX_FILENAME, header=None)
    tax.columns = ['Refseq ID', 'Taxon ID']
    ccs = ccs.merge(
        tax,
        left_on='Refseq ID',
        right_on='Refseq ID',
        how='outer'
    )
    md5sums = pd.read_csv(input_directory / MD5SUM_FILE, sep='\t', header=None).iloc[:, :2]
    md5sums.columns = ['Identifier', 'MD5 Sum']
    md5sums['Refseq ID'] = md5sums['Identifier'].apply(
        lambda x: filename_id_mapper[Path(x).name[:-3]]
    )
    md5sums.pop('Identifier')
    ccs = ccs.merge(
        md5sums,
        left_on='Refseq ID',
        right_on='Refseq ID',
        how='outer'
    )
    return preds, expl, ccs


@click.command()
@click.argument('input_directory', type=click.Path(file_okay=False, exists=True))
def preprocess(input_directory):
    preds, expl, ccs = load_data(input_directory)
    preds.reset_index(drop=True).to_feather('predictions.feather')
    expl.reset_index(drop=True).to_feather('explanations.feather')
    ccs.reset_index(drop=True).to_feather('compleconta.feather')


if __name__ == '__main__':
    preprocess()