#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 08/05/2018.
#
import os
import sys
import shutil
import ftplib
import tempfile
import argparse
from time import sleep

import django
from Bio import Entrez
from datetime import date, timedelta, datetime
from redis import Redis
from rq import Queue

from phenotypePredictionApp.variables import PHENDB_BASEDIR, PHENDB_QUEUE, PHENDB_DEBUG
from enqueue_job import phenDB_enqueue

parser = argparse.ArgumentParser()
parser.add_argument("-d", "--days_back", default=60, help="Precalculate sequences from refseq released up to days in the past")
parser.add_argument("-l", "--latest", default=None, help="Latest release date of sequence to precalculate (format: YYYY/MM/DD)")
parser.add_argument("-m", "--max_n", default=None, help="Maximum number of sequences to precalculate")
args = parser.parse_args()

Entrez.email = "lukas.lueftinger@univie.ac.at"


def get_latest_refseq_genomes(n_days, only_reference=False, max_n=None, latest=None):
    records = []
    latest_date = date.today() if latest is None else datetime.strptime(latest, "%Y/%m/%d").date()
    minus_n_days = latest_date - timedelta(days=int(n_days))
    dateformatted = minus_n_days.strftime("%Y/%m/%d")
    latestformatted = latest_date.strftime("%Y/%m/%d")
    representative_str = "" if only_reference else "OR \"representative genome\"[filter]"
    search_string = 'bacteria[filter] ' \
                    'AND ("reference genome"[filter] {rs}) ' \
                    'AND ("{minus}"[SeqReleaseDate] : "{today}"[SeqReleaseDate])'.format(rs=representative_str,
                                                                                         minus=dateformatted,
                                                                                         today=latestformatted)
    with Entrez.esearch(db="assembly", term=search_string, retmax=9999999) as handle:
        record = Entrez.read(handle)
    idlist = record["IdList"]
    print("Number of entries found: ", len(idlist))
    for i in idlist:
        retrycount = 5
        for times in range(retrycount):
            try:
                with Entrez.esummary(db="assembly", id=i) as summary_handle:
                    summary = Entrez.read(summary_handle, validate=False)
                    summarydict = summary["DocumentSummarySet"]["DocumentSummary"][-1]
                    taxid = summarydict["Taxid"]
                    name = summarydict["SpeciesName"]
                    assembly_id = summarydict["LastMajorReleaseAccession"]
                    ftppath = summarydict["FtpPath_RefSeq"]
                    if ftppath != "":
                        records.append((name, taxid, assembly_id, ftppath))
                break
            except Exception as e:
                print(e)
                print("Retrying {n} more times...".format(n=retrycount - times))
    if max_n is not None:
        records = records[:max_n]
    for entry in records:
        print("Retrieved RefSeq entry for {e}.".format(e=entry[0]))
    return records


def download_genomes(los, path):
    if os.path.exists(path):
        shutil.rmtree(path)
    with tempfile.TemporaryDirectory() as tmpname:
        for name, taxid, assembly_id, ftppath in los:
            ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov", "anonymous", "password")
            restpath = "/".join(ftppath.split("/")[3:])
            ftp.cwd("/{rp}".format(rp=restpath))
            genomicfile = [x for x in ftp.nlst() if "genomic.fna.gz" in x and "from" not in x][0]
            if not genomicfile:
                print("No genome found for {n}".format(n=name))
                continue
            fullpath_local = os.path.join(tmpname, genomicfile)
            with open(fullpath_local, "wb") as outfile:
                ftp.retrbinary("RETR {file}".format(file=genomicfile), outfile.write)
            shutil.move(fullpath_local, os.path.join(tmpname, "PHENDB_PRECALC_" + assembly_id + ".fna.gz"))
            print("Downloaded genome FASTA for {n}.".format(n=name))
        shutil.copytree(tmpname, path)


def check_add_precalc_job():
    from phenotypePredictionApp.models import Job
    try:
        Job.objects.get(key="PHENDB_PRECALC")
    except:
        new_precalc_job = Job(key="PHENDB_PRECALC")
        new_precalc_job.save()


def add_taxids_to_precalc_bins(los):
    from phenotypePredictionApp.models import Bin, Job, Taxon
    for name, taxid, assembly_id, ftppath in los:
        binname = "PHENDB_PRECALC_" + assembly_id + ".fna.gz"
        givenbin = Bin.objects.filter(bin_name=binname)
        given_taxon = Taxon.objects.filter(tax_id=taxid)
        if givenbin and given_taxon:
            givenbin.update(tax_id=str(taxid),
                            assembly_id=str(assembly_id),
                            taxon_name=str(given_taxon[0].taxon_name),
                            taxon_rank=str(given_taxon[0].taxon_rank))


def save_unadded_genome_ids(savepath, los):
    with open(savepath, "w") as outfile:
        for entry in los:
            outfile.write("\t".join(entry))
            outfile.write("\n")


def load_unadded_genome_ids(savepath):
    if os.path.exists(savepath):
        with open(savepath, "r") as infile:
            los = [x.strip().split("\t") for x in infile.readlines()]
        os.remove(savepath)
        return los
    return None


def main():

    unadded_saveloc = os.path.join(PHENDB_BASEDIR, "logs/unadded_refseq_genomes.tsv")
    ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
    infolder = os.path.join(PHENDB_BASEDIR, "data/uploads/PHENDB_PRECALC")
    outfolder = os.path.join(PHENDB_BASEDIR, "data/results/PHENDB_PRECALC_results")
    logfolder = os.path.join(outfolder, "logs")
    pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/picaPipeline.nf")

    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    os.environ["PYTHONPATH"] = ppath
    django.setup()

    os.makedirs(outfolder, exist_ok=True)
    os.makedirs(logfolder, exist_ok=True)

    print("Checking if we have unadded refseq entries sitting around...")
    unadded_los = load_unadded_genome_ids(unadded_saveloc)
    if unadded_los:
        print("Found unadded entries. Using these instead of requesting new ones.")
        gtlist = unadded_los
    else:
        print("Downloading list of genomes from RefSeq...")
        gtlist = get_latest_refseq_genomes(n_days=args.days_back, latest=args.latest, max_n=args.max_n)

    while len(gtlist) < 0:
        print("To end precalculation now, press Ctrl+C within the next 30 sec.")
        try:
            sleep(30)
        except KeyboardInterrupt:
            print("Precalculation interrupted. "
                  "{n} entries will be written to the backup file to be added later.".format(n=len(gtlist)))
            save_unadded_genome_ids(savepath=unadded_saveloc, los=gtlist)
            os._exit(0)

        fifty = gtlist[-50:]
        gtlist = gtlist[:-50]

        print("Downloading and processing <= 50 genomes...")
        print("Backlog of sequences: {n}".format(n=len(gtlist)))
        download_genomes(los=fifty, path=infolder)
        print("Submitting precalculation job. Bins in folder {inf} will be added to the database.".format(inf=infolder))
        check_add_precalc_job()
        q = Queue(PHENDB_QUEUE, connection=Redis())
        pipeline_call = q.enqueue_call(func=phenDB_enqueue,
                                       args=(ppath, pipeline_path, infolder, outfolder, 0.5, ""),
                                       timeout='72h',
                                       ttl='72h',
                                       )
        while pipeline_call.result is None:
            sleep(30)

        if pipeline_call.result is 0:
            print("Precalculation was successful.")
            print("Adding taxonomic information to precalculated bins.")
            add_taxids_to_precalc_bins(fifty)
            print("Batch finished. added {lolos} items to database.".format(lolos=len(fifty)))
        else:
            print("Precalculation has failed!")

if __name__ == "__main__":
    main()
