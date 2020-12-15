#!/usr/bin/env python3

import shutil

from Bio import SeqIO
from Bio.Alphabet import IUPAC, HasStopCodon
from Bio import Alphabet
import sys, os


full_prot_ab = Alphabet._consensus_alphabet([IUPAC.extended_protein, HasStopCodon(IUPAC.protein)])


def parse_FASTA(inputfasta, binname, errorfile):
    isdna = True
    isprotein = True

    with open("sanitychecked.fasta", "w") as outfile:
        for read in SeqIO.parse(inputfasta, "fasta", IUPAC.ambiguous_dna):
            if not Alphabet._verify_alphabet(read.seq):
                isdna = False
                break
            SeqIO.write(read, outfile, "fasta")
    if isdna and os.stat("sanitychecked.fasta").st_size != 0:
        print("DNA", end="")
        shutil.copy(inputfasta, 'sanitychecked.fasta')
        return

    with open("sanitychecked.fasta", "w") as outfile:
        for read in SeqIO.parse(inputfasta, "fasta", full_prot_ab):
            if not Alphabet._verify_alphabet(read.seq):
                isprotein = False
                break
            SeqIO.write(read, outfile, "fasta")
    if isprotein and os.stat("sanitychecked.fasta").st_size != 0:
        print("protein", end="")
        shutil.copy(inputfasta, 'sanitychecked.fasta')
        return

    if (not isdna and not isprotein) or os.stat("sanitychecked.fasta").st_size == 0:
        with open(errorfile, "a") as myfile:
            myfile.write(
                "WARNING: The file {b} is empty or not a valid fasta file "
                "and was dropped from the analysis.\n\n".format(b=binname)
            )
        os.remove("sanitychecked.fasta")
        os._exit(1)


if __name__ == '__main__':
    binname = sys.argv[1]
    item = sys.argv[2]
    errorfile = sys.argv[3]
    parse_FASTA(item, binname, errorfile)
