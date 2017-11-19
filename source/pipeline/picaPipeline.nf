#!/usr/bin/env nextflow
// define static variables

file(params.workdir).mkdir()
jobname = file(params.inputfolder).getBaseName()
outdir = "$params.workdir/${jobname}_results/"
file(outdir).mkdir()
models = file(params.modelfolder).listFiles()
input_files = Channel.fromPath("${params.inputfolder}/*.fasta")
input_gzipfiles = Channel.fromPath("${params.inputfolder}/*.tar.gz")
all_input_files = Channel.fromPath("${params.inputfolder}/*")
hmmdb = file(params.hmmdb)

log.info"""
    ##################################################################################
    
    PhenDB Pipeline started.

    Input folder containing bin files (--inputfolder): $params.inputfolder
    Output directory: $outdir
    Job name: $jobname
    
    Disabled compute nodes (for hmmer computation) (--omit_nodes): ${(params.omit_nodes == "") ? "None" : params.omit_nodes }
    Accuracy cutoff for displaying PICA results (--accuracy_cutoff): $params.accuracy_cutoff

    ##################################################################################
    """.stripIndent()

// initialize progressbar file
progressbar_hmmer = file("$outdir/logs/progress.log")
progressbar_hmmer.text = ""

// initialize filecount file
fastafilecount= file("$outdir/logs/fastafilecount.log")
fastafilecount.text = ""

// initialize error file for "sanity" errors (eg. corrupt fasta files)
errorfile= file("$outdir/logs/sanity_errors.log")
errorfile.text=""

// unzip tar.gz files
process tgz {

    input:
    file(tarfile) from input_gzipfiles

    output:
    file("${tarfile.getSimpleName()}/*.fasta") into tgz_unraveled_fasta
    file("${tarfile.getSimpleName()}/*") into tgz_unraveled_all

    script:
    """
    tar -xf $tarfile
    """
}

// combine raw fasta files and those extracted from tar.gz files
all_fasta_input_files = input_files.mix(tgz_unraveled_fasta.flatten())
all_fasta_input_files.into{all_fasta_input_files1; all_fasta_input_files2}
truly_all_input_files = all_input_files.mix(tgz_unraveled_all.flatten()) // route this through file extension checking map


// Error handling
truly_all_input_files.subscribe {

    //check if there are any non-fasta files
    if ((!(it.getName() ==~ /.+\.fasta$/ )) && (!(it.getName() =~ /.+\.tar.gz$/ ))) {

        //TODO: print to logger
        println "WARNING: Some of the files you uploaded do not end with '.fasta' or '.tar.gz'- I cannot analyze them"
        errorfile.append("WARNING: Some of the files you uploaded do not end with '.fasta' or '.tar.gz'- I cannot analyze them \n")
    }
    //check if there are any files with non-ascii character names
    if (!( it.getBaseName() ==~ /^\p{ASCII}+$/ )) {
        //TODO: print to logger
        println "WARNING: Some of the files you uploaded include non-ASCII characters - Analysis may fail!"
    }
}

// Passes every fasta file through Biopythons SeqIO to check for corrupted files    TODO: test this
process fasta_sanity_check {
    errorStrategy 'ignore'

    input:
    file(item) from all_fasta_input_files1

    output:
    set val(binname), file("sanitychecked.fasta") into fasta_sanitycheck_out

    script:
    binname = item.getBaseName()
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.Alphabet import IUPAC
    from Bio import Alphabet
    import sys, os

    
    with open("sanitychecked.fasta","w") as outfile:
        for read in SeqIO.parse("${item}", "fasta", IUPAC.ambiguous_dna):
            if not Alphabet._verify_alphabet(read.seq):
                with open("${errorfile}", "a") as myfile:
                    myfile.write("There was an unexpected letter in the sequence, aborted further processing. Allowed letters are G,A,T,C,R,Y,W,S,M,K,H,B,V,D,N ")
                os.remove("sanitychecked.fasta")
                sys.exit("There was an unexpected letter in the sequence, aborting")
            SeqIO.write(read, outfile, "fasta")

    """

}

// Print the number of fasta files to a file for progress display
all_fasta_input_files2.count().subscribe {fastafilecount.text=it} //LL: beautiful!

process md5sum {

    tag { binname }

    input:
    set val(binname), file(item) from fasta_sanitycheck_out

    output:
    set val(binname), stdout, file(item) into md5_out

    script:
    """
    echo -n \$(md5sum ${item} | cut -f1 -d" ")
    """
}

// call prodigal for every sample in parallel
// output each result as a set of the sample id and the path to the prodigal outfile
process prodigal {

    tag { binname }

    module "prodigal"
    memory = "2 GB"

    input:
    set val(binname), val(mdsum), file(item) from md5_out

    output:
    set val(binname), val(mdsum), file("prodigalout.faa") into prodigalout

    script:
    """
    prodigal -i ${item} -a prodigalout.faa > /dev/null
    """
}

// call hmmer daemon for every sample in series
// very un-nextflow, but necessary: after each completed round, the ID is written to global file progressbar_hmmer
// ((number of lines in file/number of samples)*100)) == overall progress of hmmer
process hmmer {

    tag { binname }

    maxForks 1  //do not parallelize!
    module "hmmer"

    input:
    set val(binname), val(mdsum), file(item) from prodigalout

    output:
    set val(binname), val(mdsum), file("hmmer.out"), file(item) into hmmerout

    script:
    """
    DISABLED=\$(echo ${params.omit_nodes} | sed -e 's/\\([^ ]\\+\\)/-e &/g')
    if [[ -n "\$DISABLED" ]] ; then
        HMM_DAEMONCLIENTS=\$(echo cubeb{01..30} | tr " " "\\n" | grep -v \$(echo \$DISABLED) | tr "\\n" " ")
    else
        HMM_DAEMONCLIENTS=\$(echo cubeb{01..30})
    fi
    echo \$HMM_DAEMONCLIENTS
    
    hmmc.py -i ${item} -d $hmmdb -s \$HMM_DAEMONCLIENTS -n 100 -q 5 -m 1 -o hmmer.out
    echo ${binname} >> ${progressbar_hmmer}
    """
}

// compute contamination and completeness using compleconta
process compleconta {

    tag { binname }

    module "muscle"
    module "compleconta/0.1"

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) from hmmerout

    output:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem), file("complecontaitem.txt") into complecontaout

    """
    compleconta.py $prodigalitem $hmmeritem | tail -1 > complecontaitem.txt
    """
}


complecontaout.into{complecontaout_continue; bin_to_db}

// write the bin properties to database
//TODO: disable passing of errors when adding result_enog rows after we have added all enogs to database
process write_bin_to_db { //TODO: implement checking if bin already exists

    tag { binname }

    //errorStrategy 'ignore'

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem), file(complecontaitem) from bin_to_db

    script:
"""
#!/home/user/lueftinger/miniconda3/envs/py3env/bin/python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    

django.setup()
from phenotypePredictionApp.models import *

try:
    parentjob = job.objects.get(job_name="${jobname}")
except ObjectDoesNotExist:
    sys.exit("Job not found.")

# get completeness and contamination
with open("${complecontaitem}", "r") as ccfile:
    cc = ccfile.readline().split()

# write bin to database
try:
    newbin = bin(bin_name="${binname}",
                 comple=cc[0],
                 conta=cc[1],
                 md5sum="${mdsum}",
                 job=job.objects.get(job_name="${jobname}"))
    
    newbin.save()
except IntegrityError:
    sys.exit("Skipping adding of already known bin")
    
#write hmmer output to result_enog table
with open("${hmmeritem}", "r") as enogresfile:
    for line in enogresfile:
        try:
            enogfield = line.split()[1]
            print("Attempting to add enog {en} from bin {bn} to database.".format(en=enogfield,
                                                                                  bn="${binname}"))
            newenog = result_enog(bin=bin.objects.get(bin_name="${binname}"),
                                  enog=enog.objects.get(enog_name=enogfield))
            newenog.save()
        except IntegrityError:
            print("Skipping due to IntegrityError.")
            continue               

"""
}

// compute accuracy from compleconta output and model intrinsics (once for each model).
process accuracy {

    tag { "${binname}_${model.getBaseName()}" }

    memory = '10 MB'
    errorStrategy 'ignore'  //model files not yet complete, TODO: remove this!!!!

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem), file(complecontaitem) from complecontaout_continue
    each model from models

    output:
    set val(binname), val(mdsum), val(model), file(hmmeritem), file(prodigalitem), file(complecontaitem), stdout into accuracyout

    when:
    model.isDirectory()

    script:
    RULEBOOK = model.getBaseName()
    ACCURACYFILE = "$model/${RULEBOOK}.accuracy"
    """
    python2 $params.balanced_accuracy_path $ACCURACYFILE $complecontaitem
    """
}

// call pica for every sample for every condition in parallel
process pica {

    tag { "${binname}_${model.getBaseName()}" }

    module 'pica'
    memory = '500 MB'
    errorStrategy 'ignore'  //model files not yet complete, TODO: remove this!!!!

    input:
    set val(binname), val(mdsum), val(model), file(hmmeritem), file(prodigalitem), file(complecontaitem), val(accuracy) from accuracyout

    output:
    set val(binname), val(mdsum), val(RULEBOOK), stdout, val(accuracy) into picaout  //print decision on stdout, and put stdout into return set

    script:
    RULEBOOK = model.getBaseName()
    TEST_MODEL = "$model/${RULEBOOK}.rules"
    float accuracy_cutoff = params.accuracy_cutoff as float
    float accuracy_float = accuracy as float

    if (accuracy_float >= accuracy_cutoff) {
    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f2 $hmmeritem | tr "\n" "\t" >> tempfile.tmp
    test.py -m $TEST_MODEL -t $RULEBOOK -s tempfile.tmp > picaout.result
    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3 | tr "\t" " ")
    """
    }

    else {
    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f2 $hmmeritem | tr "\n" "\t" >> tempfile.tmp
    echo "N/A\tN/A" > picaout.result
    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3 | tr "\t" " ")
    """
    }
}

// merge all results into a file called $id.results and move each file to results folder.
picaout.into{pica_db_write; pica_out_write}

outfilechannel = pica_out_write.collectFile() { item ->
    [ "${item[0]}.results", "${item[2]} ${item[3]} ${item[4]}" ]  // use given bin name as filename
}

process tar_results {

    tag { jobname }

    stageInMode 'copy'  //actually copy in the results so we not only tar symlinks

    input:
    val(allfiles) from outfilechannel.collect()

    output:
    file("${jobname}.tar.gz") into tgz_out

    script:
    """
    mkdir ${jobname}
    mv ./*.results ${jobname}
    tar -cvf ${jobname}.tar.gz ./${jobname}
    """
}

tgz_out.into { tgz_write_out; tgz_to_db }
tgz_write_out.subscribe { it.copyTo(outdir) }

db_written = pica_db_write.collectFile() { item ->
    [ "${item[1]}.results", "${item[2]} ${item[3]} ${item[4]}" ]  // use md5sum as filename
}

process write_pica_result_to_db { //TODO: change python executable when migrating to vm!

    //errorStrategy 'ignore'

    input:
    file(mdsum_file) from db_written

    output:
    stdout exo

    script:
"""
#!/usr/bin/env python

###n## !/home/user/lueftinger/miniconda3/envs/py3env/bin/python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import job, bin, model, result_model

# get job from db
try:
    parentjob = job.objects.get(job_name="${jobname}")
except ObjectDoesNotExist:
    sys.exit("Job not found.")

# get bin from db
try:
    parentbin = bin.objects.get(md5sum="${mdsum_file.getBaseName()}")
except ObjectDoesNotExist:
    sys.exit("Bin for this result not found.")

conditions = []
with open("${mdsum_file}", "r") as picaresults:
    for line in picaresults:
        conditions.append(line.split())

for result in conditions:
    try:
        get_bool = {"YES": True, "NO": False, "N/A": None}
        boolean_verdict = get_bool[result[1]]
        #get model from db
        try:
            this_model=model.objects.get(model_name=result[0], is_newest=True)
        except ObjectDoesNotExist:
            sys.exit("Current Model for this result not found.")
        
        modelresult = result_model(verdict=boolean_verdict,
                                   accuracy=float(result[3].split()[0]),
                                   bin=parentbin,
                                   model=this_model
                                   )
        modelresult.save()
    except (IntegrityError, ) as e:
        sys.exit(e)
"""
}

process write_tgz_to_db {

    input:
    file(tgz) from tgz_to_db

    script:
"""
#!/home/user/lueftinger/miniconda3/envs/py3env/bin/python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    

django.setup()
from phenotypePredictionApp.models import ResultFile

try:
    new_tgz = ResultFile(actualID=${jobname},
                         document=${tgz}
                         )
except IntegrityError:
    sys.exit("Exited with integrity error upon adding results to database.")
"""
}






workflow.onComplete {
    println "picaPipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
