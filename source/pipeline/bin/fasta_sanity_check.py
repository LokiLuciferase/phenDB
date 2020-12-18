#!/usr/bin/env python3
from typing import List
import shutil

from Bio import SeqIO, Seq
from Bio.Alphabet import IUPAC, HasStopCodon
from Bio import Alphabet
import sys, os


MIN_TOTAL_DNA_LEN = 20_000
full_prot_ab = Alphabet._consensus_alphabet([IUPAC.extended_protein, HasStopCodon(IUPAC.protein)])


def try_read_fasta(inputfile: str) -> List[str]:
    seqs = list(SeqIO.parse(inputfile, 'fasta'))
    if not any(seqs):
        raise ValueError('No sequences found.')
    return [str(x.seq) for x in seqs]


def is_valid_dna_seqs(seqs: List[str]) -> bool:
    return all(Alphabet._verify_alphabet(Seq.Seq(x, IUPAC.ambiguous_dna)) for x in seqs)


def is_valid_protein_seqs(seqs: List[str]) -> bool:
    return all(Alphabet._verify_alphabet(Seq.Seq(x, full_prot_ab)) for x in seqs)


def dna_size_ok(seqs: List[str]) -> bool:
    return (
            sum(len(x) for x in seqs) > MIN_TOTAL_DNA_LEN
            and all(len(x) > 0 for x in seqs)
    )


def fail_with_log(error: str, error_file: str):
    with open(error_file, 'a') as fout:
        fout.write(error + '\n\n')
    os._exit(1)
    raise RuntimeError(error)


def main(input_file: str, bin_name: str, error_file: str):
    not_valid_seq_err = f'WARNING: The file {bin_name} is empty or not a valid fasta file ' \
                        f'and was dropped from the analysis.'
    not_long_enough_err = f'WARNING: The file {bin_name} was too small ' \
                          f'and was dropped from the analysis.'
    try:
        seqs = try_read_fasta(input_file)
    except:
        seqs = None
        fail_with_log(not_valid_seq_err, error_file)

    file_type = None
    if is_valid_dna_seqs(seqs):
        if dna_size_ok(seqs):
            file_type = 'DNA'
        else:
            fail_with_log(not_long_enough_err, error_file)
    elif is_valid_protein_seqs(seqs):
        file_type = 'protein'
    else:
        pass

    if file_type is None:
        fail_with_log(not_valid_seq_err, error_file)
    else:
        print(file_type, end='')
        shutil.copy(input_file, 'sanitychecked.fasta')


if __name__ == '__main__':
    main(input_file=sys.argv[1], bin_name=sys.argv[2], error_file=sys.argv[3])
