#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 08/05/2018.
#
import os
import sys
import glob
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

from phenotypePredictionApp.variables import (
    PHENDB_BASEDIR,
    PHENDB_DATA_DIR,
    PHENDB_QUEUE,
    PHENDB_DEBUG,
)
from enqueue_job import phenDB_enqueue, phenDB_recalc

django.setup()
from phenotypePredictionApp.models import Bin, Job, BinInJob, Taxon


parser = argparse.ArgumentParser()
parser.add_argument(
    "-d",
    "--days_back",
    default=60,
    help="Precalculate sequences " "from refseq released up to days in the past",
)
parser.add_argument(
    "-l",
    "--latest",
    default=None,
    help="Latest release date of sequence " "to precalculate (format: YYYY/MM/DD)",
)
parser.add_argument(
    "-b", "--n_batches", default=999999, help="How many 50-genome batches to process until stopping"
)
parser.add_argument(
    "-m", "--max_n", default=None, help="Maximum number of sequences to precalculate"
)
parser.add_argument(
    "-t",
    "--update_taxids",
    default=False,
    help="Only get list of refseq entries, " "and add taxonomy info to those present in the DB.",
)
parser.add_argument(
    "-r",
    "--rerun_existing",
    default=False,
    action="store_true",
    help="Re-enter known genomes to PhenDB (for re-calculation of results with updated models)",
)
parser.add_argument(
    "-e",
    "--get_explanations",
    default=False,
    action="store_true",
    help="Whether to compute SHAP explanations for precalculated bins."
)
args = parser.parse_args()

Entrez.email = "lukas.lueftinger@univie.ac.at"
REFSEQ_GENOMES_BACKUP_LOC =  f"{PHENDB_DATA_DIR}/refseq_genomes_precalc"
RECALC_MAX_BATCH_NO = 100


# download accession IDs, names, taxon IDs and FTP paths for refseq genomes submitted in the given time span
def get_latest_refseq_genomes(n_days, only_reference=False, max_n=None, latest=None):
    records = []
    latest_date = date.today() if latest is None else datetime.strptime(latest, "%Y/%m/%d").date()
    minus_n_days = latest_date - timedelta(days=int(n_days))
    dateformatted = minus_n_days.strftime("%Y/%m/%d")
    latestformatted = latest_date.strftime("%Y/%m/%d")
    representative_str = "" if only_reference else 'OR "representative genome"[filter]'
    search_string = (
        "bacteria[filter] "
        'AND ("reference genome"[filter] {rs}) '
        'AND ("{minus}"[SeqReleaseDate] : "{today}"[SeqReleaseDate])'.format(
            rs=representative_str, minus=dateformatted, today=latestformatted
        )
    )
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


# download refseq genomes given in list los and save them to folder path
def download_genomes(los, path):
    if os.path.exists(path):
        for genome in glob.glob(os.path.join(path, "*.fna.gz")):
            try:
                shutil.move(genome, REFSEQ_GENOMES_BACKUP_LOC)
            except shutil.Error:
                pass
        shutil.rmtree(path)
    with tempfile.TemporaryDirectory() as tmpname:
        for name, taxid, assembly_id, ftppath in los:
            finished = False
            while not finished:
                try:
                    ftp = ftplib.FTP("ftp.ncbi.nlm.nih.gov", "anonymous", "password")
                    restpath = "/".join(ftppath.split("/")[3:])
                    ftp.cwd("/{rp}".format(rp=restpath))
                    genomicfile = [
                        x for x in ftp.nlst() if "genomic.fna.gz" in x and "from" not in x
                    ][0]
                    if not genomicfile:
                        print("No genome found for {n}".format(n=name))
                        continue
                    fullpath_local = os.path.join(tmpname, genomicfile)
                    with open(fullpath_local, "wb") as outfile:
                        ftp.retrbinary("RETR {file}".format(file=genomicfile), outfile.write)
                    shutil.move(
                        fullpath_local,
                        os.path.join(tmpname, "PHENDB_PRECALC_" + assembly_id + ".fna.gz"),
                    )
                    print("Downloaded genome FASTA for {n}.".format(n=name))
                    finished = True
                except ftplib.error_temp:
                    print("Trying again after FTP temporary error. Sleeping 30 sec first.")
                    sleep(5)
                except EOFError:
                    print("Trying again after FTP temporary error. Sleeping 30 sec first.")
                    sleep(5)
        shutil.copytree(tmpname, path)


# check if a job with id "PHENDB_PRECALC" exists in the DB; if not, create it (required for phenDB pipeline)
def check_add_precalc_job(jobname="PHENDB_PRECALC"):
    print("Checking if PHENDB_PRECALC job exists in database...")
    try:
        pc_job = Job.objects.get(key=jobname)
        pc_job.job_date = datetime.now()
        pc_job.save()
        print("PHENDB_PRECALC job date updated.")
    except:
        new_precalc_job = Job(key=jobname)
        new_precalc_job.save()
        print("PHENDB_PRECALC job created.")


# attempt to add scientific names, taxids and accession ids to bins with the given name in the DB. Using Taxon table
def add_taxids_to_precalc_bins(los):
    precalc_job = Job.objects.get(key="PHENDB_PRECALC")
    for name, taxid, assembly_id, ftppath in los:
        bin_alias = "PHENDB_PRECALC_" + assembly_id + ".fna.gz"
        givenbij = BinInJob.objects.filter(bin_alias=bin_alias, job=precalc_job)
        given_taxon = Taxon.objects.filter(tax_id=taxid)
        if givenbij and given_taxon:
            givenbin = givenbij[0].bin
            givenbin.tax_id = str(taxid)
            givenbin.assembly_id = str(assembly_id)
            givenbin.save()


# if precalculation is interrupted at save point, save all remaining refseq FTP paths and ids to a file
def save_unadded_genome_ids(savepath, los):
    with open(savepath, "w") as outfile:
        for entry in los:
            outfile.write("\t".join(entry))
            outfile.write("\n")


# check if a savefile exists for unadded refseq paths, if so, load it and use them to download their genomes
def load_unadded_genome_ids(savepath):
    if os.path.exists(savepath):
        with open(savepath, "r") as infile:
            los = [x.strip().split("\t") for x in infile.readlines()]
        os.remove(savepath)
        return los
    return None


# Re-add genomes from local cache to phenDB to make predictions on new phenotrex models
def rerun_known_genomes(ppath, outfolder):
    for batch_no in range(RECALC_MAX_BATCH_NO):
        print("Batch {i}: Starting...".format(i=batch_no))
        check_add_precalc_job()
        pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/main.nf")
        q = Queue(PHENDB_QUEUE, connection=Redis())
        pipeline_call = q.enqueue_call(
            func=phenDB_recalc,
            args=(ppath, pipeline_path, outfolder, batch_no, RECALC_MAX_BATCH_NO),
            timeout="72h",
            ttl="72h",
        )
        while pipeline_call.result is None:
            sleep(30)

        if pipeline_call.result == 0:
            print("Batch {i}: Recalculation was successful.".format(i=batch_no))
    print("All recalculations successful.")
    return True


# submit a PhenDB job to redis queue, with job ID "PHENDB_PRECALC", sleep until finished, return True if success
def start_precalc_queue(ppath, infolder, outfolder):
    os.environ['PYTHONPATH'] = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
    from enqueue_job import phenDB_enqueue, phenDB_recalc

    pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/main.nf")
    q = Queue(PHENDB_QUEUE, connection=Redis())
    pipeline_call = q.enqueue_call(
        func=phenDB_enqueue,
        args=(ppath, pipeline_path, infolder, outfolder, "", args.get_explanations),
        timeout="72h",
        ttl="72h",
    )
    while pipeline_call.result is None:
        sleep(30)

    if pipeline_call.result == 0:
        print("Precalculation was successful.")
        return True
    print("Precalculation has failed!")
    return False


# create in, out and logfolders for phenDB input, download desired refseq genomes,
# split into chunks to be processed sequentially
def main():
    unadded_saveloc = os.path.join(PHENDB_DATA_DIR, "logs/unadded_refseq_genomes.tsv")
    infolder = os.path.join(PHENDB_DATA_DIR, "uploads/PHENDB_PRECALC")
    outfolder = os.path.join(PHENDB_DATA_DIR, "results/PHENDB_PRECALC_results")
    ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
    logfolder = os.path.join(outfolder, "logs")

    os.makedirs(outfolder, exist_ok=True)
    os.makedirs(logfolder, exist_ok=True)

    # handle special precalculation cases
    if args.update_taxids:
        print("Only adding taxonomic info to existing bins, then exiting.")
        gtlist = get_latest_refseq_genomes(n_days=args.days_back, latest=args.latest)
        add_taxids_to_precalc_bins(los=gtlist)
        os._exit(0)
    elif args.rerun_existing:
        print("Re-running old refseq genomes with new models.")
        rerun_known_genomes(ppath=ppath, outfolder=outfolder)
        os._exit(0)

    # handle default case
    print("Checking if we have unadded refseq entries sitting around...")
    unadded_los = load_unadded_genome_ids(unadded_saveloc)
    if unadded_los:
        print("Found unadded entries. Using these instead of requesting new ones.")
        gtlist = unadded_los
    else:
        print("Downloading list of genomes from RefSeq...")
        gtlist = get_latest_refseq_genomes(
            n_days=args.days_back, latest=args.latest, max_n=args.max_n
        )

    init_batches = int(args.n_batches)
    while len(gtlist) > 0:
        save_unadded_genome_ids(savepath=unadded_saveloc, los=gtlist)
        if init_batches <= 0:
            sys.stdout.write(
                "Finishing after reaching max number of batches to process.\n"
                "{n} entries will be written to the backup file to be added later.\n".format(
                    n=len(gtlist)
                )
            )
            sys.stdout.flush()
            os._exit(0)

        sys.stdout.write(
            "Remaining genome IDs saved. To end precalculation now, "
            "press Ctrl+C within the next 30 sec.\n"
        )
        sys.stdout.flush()
        try:
            sleep(30)
        except KeyboardInterrupt:
            os._exit(0)

        fifty = gtlist[-50:]
        gtlist = gtlist[:-50]

        print("Downloading and processing <= 50 genomes...")
        print("Backlog of sequence IDs: {n}".format(n=len(gtlist)))
        download_genomes(los=fifty, path=infolder)
        print(
            "Submitting precalculation job. Bins in folder {inf} will be added to the database.".format(
                inf=infolder
            )
        )
        check_add_precalc_job()
        exitstat = start_precalc_queue(ppath=ppath, infolder=infolder, outfolder=outfolder)
        if exitstat:
            print("Adding taxonomic information to precalculated bins.")
            add_taxids_to_precalc_bins(fifty)
            print("Batch finished. added {lolos} items to database.".format(lolos=len(fifty)))
            init_batches -= 1


if __name__ == "__main__":
    main()
