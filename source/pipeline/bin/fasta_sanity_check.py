#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Alphabet import IUPAC, HasStopCodon
from Bio import Alphabet
import sys, os
import gzip

binname = sys.argv[1]
item = sys.argv[2]
errorfile = sys.argv[3]

full_prot_ab = Alphabet._consensus_alphabet([IUPAC.extended_protein, HasStopCodon(IUPAC.protein)])


def parse_FASTA(inputfasta):
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
        return

    with open("sanitychecked.fasta", "w") as outfile:
        for read in SeqIO.parse(inputfasta, "fasta", full_prot_ab):
            if not Alphabet._verify_alphabet(read.seq):
                isprotein = False
                break
            SeqIO.write(read, outfile, "fasta")
    if isprotein and os.stat("sanitychecked.fasta").st_size != 0:
        print("protein", end="")
        return

    if (not isdna and not isprotein) or os.stat("sanitychecked.fasta").st_size == 0:
        with open(errorfile, "a") as myfile:
            myfile.write(
                "WARNING: The file {b} is empty or not a valid fasta file "
                "and was dropped from the analysis.\n\n".format(b=binname)
            )
        os.remove("sanitychecked.fasta")
        os._exit(1)


def invalid_gz_error():
    with open(errorfile, "a") as myfile:
        myfile.write(
            "WARNING: An error occured during parsing of {b}. "
            "The file was dropped from the analysis.\\n\\n".format(b=binname)
        )
    os.remove("sanitychecked.fasta")
    os._exit(1)


if item.endswith(".gz"):
    try:
        with gzip.open(item, "rt") as inputfasta:
            try:
                parse_FASTA(inputfasta)
            except:
                invalid_gz_error()
    except OSError:
        with open(errorfile, "a") as myfile:
            myfile.write(
                "WARNING: {b} is a compressed directory, not a file, or otherwise not a valid .gz file. "
                "It is dropped from further analysis; Nested directories cannot be analyzed! \\n\\n".format(
                    b=binname
                )
            )
            os._exit(1)
else:
    inputfasta = item
    try:
        parse_FASTA(inputfasta)
    except:
        invalid_gz_error()
