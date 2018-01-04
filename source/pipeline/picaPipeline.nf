#!/usr/bin/env nextflow
// define static variables

jobname = file(params.inputfolder).getBaseName()
outdir = file(params.outdir)
file(outdir).mkdir()
models = file(params.modelfolder).listFiles()
// make this independent of ending!! try to unpack, then directly check contents
input_files = Channel.fromPath("${params.inputfolder}/*.fasta")
input_gzipfiles = Channel.fromPath("${params.inputfolder}/*.tar.gz")
input_barezipfiles = Channel.fromPath("${params.inputfolder}/*.zip")
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
    //file("${tarfile.getSimpleName()}/*.fasta") into tgz_unraveled_fasta
    file("${tarfile.getSimpleName()}/*") into tgz_unraveled_all

    script:
    outfolder = tarfile.getSimpleName()
    """
    tar -xf $tarfile 
    rm -f $tarfile
    
    if [[ \$(ls) != "$outfolder" ]]; then
        mv ./*/ $outfolder
    fi
    """
}

// unzip .zip files
process unzip {

    input:
    file(zipfile) from input_barezipfiles

    output:
    //file("${zipfile.getSimpleName()}/*/*.fasta") into zip_unraveled_fasta
    //todo: PP: I changed this from "/*/*" to "/*" - is that ok?
    file("${zipfile.getSimpleName()}/*") into zip_unraveled_all

    script:
    outfolder = zipfile.getSimpleName()
    """
    mkdir ${outfolder} && unzip ${zipfile} -d ${outfolder}
    """
}

// combine raw fasta files and those extracted from archive files
//all_fasta_input_files = input_files.mix(tgz_unraveled_fasta.flatten(), zip_unraveled_fasta.flatten())
truly_all_input_files = all_input_files.mix(tgz_unraveled_all.flatten(), zip_unraveled_all.flatten())


//// Error handling
//truly_all_input_files.subscribe {
//
//    //check if there are any non-fasta files
//    if ((!(it.getName() ==~ /.+\.fasta$/ )) && (!(it.getName() =~ /.+\.tar.gz$/ )) && (!(it.getName() =~ /.+\.zip$/ ))){
//
//        endingmess = "WARNING: the file ${it.getName()} does not end in '.fasta', '.zip' or '.tar.gz'.\n" +
//                     "The file was dropped from the analysis.\n\n"
//        log.info(endingmess)
//        errorfile.append(endingmess)
//
//    }
//    //check if there are any files with non-ascii character names
//    if (!( it.getBaseName() ==~ /^\p{ASCII}+$/ )) {
//
//        asciimess = "WARNING: The filename of ${it.getName()} contains non-ASCII characters.\n" +
//                    "The file was dropped from the analysis.\n\n"
//        log.info(asciimess)
//        errorfile.append(asciimess)
//    }
//}

// Passes every fasta file through Biopythons SeqIO to check for corrupted files
process fasta_sanity_check {
//todo: output also sequences-only for md5sum?
    //todo: at the moment, the zip file itself is also processed here. ok for now but not pretty
    //todo: maybe recursively uncompress everything we can find and put it into one folder beforehand?
    errorStrategy 'ignore'
    tag { binname }
    maxForks 10

    input:
    file(item) from truly_all_input_files

    output:
    set val(binname), file("sanitychecked.fasta") into fasta_sanitycheck_out

    script:
    binname = item.getName()
    if (!( item.getBaseName() ==~ /^\p{ASCII}+$/ )) {

    asciimess = "WARNING: The filename of ${item.getName()} contains non-ASCII characters.\n" +
                "The file was dropped from the analysis.\n\n"
    log.info(asciimess)
    errorfile.append(asciimess)
        return //todo: this seems to raise error code 127, which does not make sense. However, it seems to serve the purpose for now
        log.info("you should not see this")
    }
// language=Python
"""
#!/usr/bin/env python
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
from Bio import Alphabet
import sys, os
import tarfile

inputfasta=""
if ("${item}".endswith("tar.gz")):
    os.makedirs("extr")
    tar = tarfile.open("${item}", "r:gz")
    for tarinfo in tar:
        #must not be another directory but a regular file
        if tarinfo.isreg():
            tar.extractall(path="extr")
            tar.close()
            inputfasta="extr/"+str(os.listdir("extr")[0])
            
        else: 
            raise 

elif ("${item}".endswith("tar")):
    os.makedirs("extr")
    tar = tarfile.open("${item}", "r:")
    for tarinfo in tar:
        if tarinfo.isreg():
            tar.extractall(path="extr")
            tar.close()
            inputfasta="extr/"+str(os.listdir("extr")[0])
    
        else: 
            raise 
else:
    inputfasta="${item}"

        

    
        
with open("sanitychecked.fasta","w") as outfile:
    for read in SeqIO.parse(inputfasta, "fasta", IUPAC.ambiguous_dna):
        if not Alphabet._verify_alphabet(read.seq):
            with open("${errorfile}", "a") as myfile:
                myfile.write("WARNING: There was an unexpected DNA letter in the sequence of file ${binname}.\\n")
                myfile.write("Allowed letters are G,A,T,C,R,Y,W,S,M,K,H,B,V,D,N.\\n")
                myfile.write("The file was dropped from the analysis.\\n\\n")
            os.remove("sanitychecked.fasta")
            sys.exit("There was an unexpected letter in the sequence, aborting.")
            
        SeqIO.write(read, outfile, "fasta")
if os.stat("sanitychecked.fasta").st_size == 0:
    with open("${errorfile}", "a") as myfile:
        myfile.write("WARNING: The file ${binname} is empty or not a fasta file and was dropped from the analysis.\\n\\n")
    os.remove("sanitychecked.fasta")
    sys.exit("The file is empty or not in fasta format, aborting.")

"""

}

// Print the number of valid fasta files to a file for progress display
fasta_sanitycheck_out.into { sanity_check_for_continue; sanity_check_for_count }
sanity_check_for_count.count().subscribe { fastafilecount.text=it }

process md5sum {

    tag { binname }

    input:
    set val(binname), file(item) from sanity_check_for_continue

    output:
    set val(binname), stdout, file(item) into md5_out

    script:
    """
    echo -n \$(md5sum ${item} | cut -f1 -d" ")
    """
}

process add_bin_to_db {

    input:
    set val(binname), val(mdsum), file(item) from md5_out

    output:
    set val(binname), val(mdsum), file(item), stdout into bin_to_db_out


    script:
// language=Python
"""
#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    

django.setup()
from phenotypePredictionApp.models import *

try:
    parentjob = UploadedFile.objects.get(key="${jobname}")
except ObjectDoesNotExist:
    sys.exit("Job not found.")

try:
    thisbin=bin.objects.get(md5sum="${mdsum}")
    print("NO", end='')
except ObjectDoesNotExist:
    # write impossible Nr. for comple and conta that would cause an error during get_accuracy if not overwritten:
    thisbin = bin(bin_name="${binname}", md5sum="${mdsum}", comple=2, conta=2)
    thisbin.save()
    print("YES", end='')

try:
    assoc= bins_in_UploadedFile(bin=thisbin, UploadedFile=parentjob)
    assoc.save()    
except IntegrityError:
    sys.exit("Cannot add bin to db: An identical file (same md5sum) from the same job is already in the db. Please remove duplicate files from your input!")
"""

}

// Do the further checks depending on if the bin is already in the db
bin_to_db_out.into{bin_is_in_db;bin_is_not_in_db}


// call prodigal for every sample in parallel
// output each result as a set of the sample id and the path to the prodigal outfile
// only execute prodigal if the been has not already been in our db
process prodigal {

    tag { binname }
    maxForks 10

    memory = "2 GB"

    input:
    set val(binname), val(mdsum), file(item), val(calc_bin_or_not) from bin_is_not_in_db

    output:
    set val(binname), val(mdsum), file("prodigalout.faa") into prodigalout

    when:
    calc_bin_or_not.equals("YES")

    script:
    """
    prodigal -i ${item} -a prodigalout.faa > /dev/null
    """
}

// if the bin has already been in our db:
process determine_models_that_need_recalculation {

    input:
    set val(binname), val(mdsum), file(item), val(calc_bin_or_not) from bin_is_in_db
    each model from models

    output:
    set val(binname), val(mdsum), val(model), stdout into determined_if_model_recalculation_needed

    when:
    calc_bin_or_not == "NO"

    script:
    // language=Python
    """
    #!/usr/bin/env python3
    
    import django
    import sys
    import os
    from django.core.exceptions import ObjectDoesNotExist
    from django.db import IntegrityError
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    
    
    django.setup()
    from phenotypePredictionApp.models import *
    
  

    try: #if this succeeds, then there is a result for this bin and the newest model version in our db
        result_model.objects.get(model=model.objects.filter(model_name="${model.getBaseName()}").latest('model_train_date'), bin=bin.objects.get(md5sum="${mdsum}"))
        print("NO", end='')    
    except:
        print("YES", end='')

    """
}

determined_if_model_recalculation_needed.into{calc_model; dont_calc_model}

// if the results for this bin and model are up-to-date, we only need to include them in the final output file
process uptodate_model_to_targz1 {

    input:
    set val(binname), val(mdsum), val(model), val(calc_model_or_not) from dont_calc_model

    output:
    set val(binname), val(mdsum), val(RULEBOOK), val(accuracy), file("verdict_and_accuracy.txt") into new_model_to_targz1_out
    // print YES or NO
    when:
    calc_model_or_not == "NO"

    script:
    RULEBOOK=model.getBaseName()
    // language=Python
    """
    #!/usr/bin/env python3
    
    import django
    import sys
    import os
    from django.core.exceptions import ObjectDoesNotExist
    from django.db import IntegrityError
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    
    
    django.setup()
    from phenotypePredictionApp.models import *
    
    #with open("verdict_and_accuracy.txt", "w") as v_a:    

    #check the database for the entries and print them to file

    """
}

process uptodate_model_to_targz2 {

    input:
    set val(binname), val(mdsum), val(RULEBOOK), val(accuracy), file(verdict_and_accuracy) from new_model_to_targz1_out
    output:
    set val(binname), val(mdsum), val(RULEBOOK), val(verdict), val(accuracy) into picaout_from_new_model //

    script:
    //TODO: read verdict and accuracy from the file
    verdict=""
    accuracy=0

}

// if the results in our db for this bin and model are outdated, a hmmerfile & compleconta file needs to be reconstructed and pica needs to be called
process old_model_to_accuracy {

    input:
    set val(binname), val(mdsum), val(model), val(calc_model_or_not) from calc_model

    output:
    set val(binname), val(mdsum), val(model), file("reconstructed_hmmer_file.txt"), file("reconstructed_compleconta_file.txt") into accuracy_in_from_old_model

    when:
    calc_model_or_not == "YES"

    script:
    // language=Python
    """
    #!/usr/bin/env python3
    
    import django
    import sys
    import os
    from django.core.exceptions import ObjectDoesNotExist
    from django.db import IntegrityError
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    
    
    django.setup()
    from phenotypePredictionApp.models import *
    
    #with open("reconstructed_hmmer_file.txt", "w") as hmmer:
        # TODO: construct hmmer file
    
    #with open("reconstructed_compleconta_file.txt", "w") as hmmer:
        # TODO: construct compleconta file
        



    """
}

// call hmmer daemon for all bins for which prodigal was run (i.e. bins that have not yet been in the db)
process hmmer {

    tag { binname }

    maxForks 1  //do not parallelize!

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
    """
}

//increase percentage completed after each hmmer job completion
// TODO: this does not work for jobs with >> 100 bins at the moment (errors in rounding)
process update_job_completeness {

    tag { jobname }

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) from hmmerout
    output:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) into job_updated_out
    script:
// language=Python
"""
#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    

django.setup()
from phenotypePredictionApp.models import *

try:
    parentjob = UploadedFile.objects.get(key="${jobname}")

    current_status = int(parentjob.job_status) if parentjob.job_status else 0
    total_valid_count = int('${fastafilecount.text}'.strip())
    plusone = int((1 / total_valid_count)*100)
    
    if (current_status + plusone < 100):
        current_status += plusone
        parentjob.job_status = current_status
        parentjob.save()
    
except ObjectDoesNotExist:
    sys.exit("Job not found.")
"""

}

// compute contamination and completeness using compleconta
process compleconta {

    tag { binname }

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) from job_updated_out

    output:
    set val(binname), val(mdsum), file(hmmeritem), file("complecontaitem.txt") into complecontaout

    """
    compleconta.py $prodigalitem $hmmeritem | tail -1 > complecontaitem.txt
    """
}


complecontaout.into{complecontaout_for_call_accuracy; complecontaout_for_hmmerresults_to_db}

process write_hmmer_results_to_db {

    tag { binname }

               // TODO: this appears to be a huge bottleneck, let's optimize this
                // TODO: check if enog.objects.in_bulk() makes sense also here

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(complecontaitem) from complecontaout_for_hmmerresults_to_db

    script:
// language=Python
"""
#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"    

django.setup()
from phenotypePredictionApp.models import *

try:
    parentbin = bin.objects.get(md5sum="${mdsum}")
except ObjectDoesNotExist:
    sys.exit("Bin not found.")


  
#write hmmer output to result_enog table
with open("${hmmeritem}", "r") as enogresfile:
    enoglist=[]
    for line in enogresfile:
        enogfield = line.split()[1]
        enoglist.append(enogfield)
    # hmmer outputs some ENOGS multiple times --> add them only once
    enoglist=list(set(enoglist))
    enog_res_objectlist=[]
    for enog_elem in enoglist:    
        try:
            print("Attempting to add enog {en} from bin {bn} to database.".format(en=enog_elem,
                                                                                  bn="${binname}"))
            new_enog_res = result_enog(bin=parentbin, enog=enog.objects.get(enog_name=enog_elem))
            enog_res_objectlist.append(new_enog_res)
        except IntegrityError:
            print("Skipping due to IntegrityError.")
            continue     
            

    try:
        result_enog.objects.bulk_create(enog_res_objectlist)        
    except IntegrityError:
        print("Could not add enogs of bin ${binname} to the db.")  
    
"""
}

process call_accuracy_for_all_models {

    tag { "${binname}_${model.getBaseName()}" }

    memory = '10 MB'

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(complecontaitem) from complecontaout_for_call_accuracy
    //todo: do not rely on folder but on db here
    each model from models

    output:
    set val(binname), val(mdsum), val(model), file(hmmeritem), file(complecontaitem) into accuracy_in

    when:
    model.isDirectory()

    script:
    """
    """
}

// compute accuracy from compleconta output and model intrinsics .
process accuracy {

    tag { "${binname}_${model.getBaseName()}" }

    memory = '10 MB'

    input:                       // .mix: this is where models from bins that were already in the db, but that have outdated results, are added
    set val(binname), val(mdsum), val(model), file(hmmeritem), file(complecontaitem) from accuracy_in.mix(accuracy_in_from_old_model)

    output:
    set val(binname), val(mdsum), val(model), file(hmmeritem), stdout into accuracyout


    script:
    """

    #!/usr/bin/env python3



    import django
    import sys
    import os
    from django.core.exceptions import ObjectDoesNotExist
    from django.db import IntegrityError
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
    django.setup()
    from phenotypePredictionApp.models import *
    import math
    def round_nearest(x, a):
        return round(round(x / a) * a, -int(math.floor(math.log10(a))))

    # get completeness and contamination
    with open("${complecontaitem}", "r") as ccfile:
        cc = ccfile.readline().split()

    try:
        parentbin = bin.objects.get(md5sum="${mdsum}")
        parentbin.comple = cc[0]
        parentbin.conta= cc[1]
        parentbin.save()


        
    except ObjectDoesNotExist:
        sys.exit("Bin not found.")
    
    #check the balanced accuracy. other metrices would be very similar to evaluate, just change the part after the last dot
    print(model_accuracies.objects.get(model=model.objects.filter(model_name="${model.getBaseName()}").latest('model_train_date'),
    comple=round_nearest(float(parentbin.comple),0.05), 
    conta=round_nearest(float(parentbin.conta),0.05)).mean_balanced_accuracy)

    """
}

// call pica for every sample for every condition in parallel
//todo: do not show YES/NO for user download-file if below threshold
process pica {

    tag { "${binname}_${model.getBaseName()}" }

    memory = '500 MB'

    input:
    set val(binname), val(mdsum), val(model), file(hmmeritem), val(accuracy) from accuracyout

    output:
    set val(binname), val(mdsum), val(RULEBOOK), stdout, val(accuracy) into picaout  //print decision on stdout, and put stdout into return set

    script:
    RULEBOOK = model.getBaseName()
    TEST_MODEL = "$model/${RULEBOOK}.rules"
//    float accuracy_cutoff = params.accuracy_cutoff as float
//    float accuracy_float = accuracy as float
//
//    if (accuracy_float >= accuracy_cutoff) {
    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f2 $hmmeritem | tr "\n" "\t" >> tempfile.tmp
    test.py -m $TEST_MODEL -t $RULEBOOK -s tempfile.tmp > picaout.result
    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3)
    """
//    }

//    else {
//    """
//    echo -ne "${binname}\t" > tempfile.tmp
//    cut -f2 $hmmeritem | tr "\n" "\t" >> tempfile.tmp
//    echo "N/A\tN/A" > picaout.result
//    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3 | tr "\t" " ")
//    """
//    }
}

// merge all results into a file called $id.results and move each file to results folder.
picaout.into{picaout_db_write; picaout_for_download}

//process replace_with_NA {
//    input:
//    set val(binname), val(mdsum), val(RULEBOOK), val(verdict), val(accuracy) from picaout_for_download.mix(picaout_from_new_model)
//
//    echo = true
//
//    output:
//    set val(binname), val(mdsum), val(RULEBOOK), stdout, val(accuracy) into NA_replaced_for_download
//
//    script:
//    float accuracy_cutoff = params.accuracy_cutoff as float
//    float accuracy_float = accuracy as float
//
//    if (accuracy_float+100 < accuracy_cutoff) {
//        newout="N/A\t"
//        newout+=verdict.split(" ")[1].toString()
//        println newout
//    }
//
//    else {
//      println verdict
//
//      }
//}

outfilechannel = picaout_for_download.mix(picaout_from_new_model).collectFile() { item ->
    [ "${item[0]}.results", "${item[2]}\t${item[3]}\t${item[4]}" ]  // use given bin name as filename
}.collect()


// create a matrix file containing all verdicts for each bin and add to output files
// create a matrix file containing summaries per model
// here we could for example add a header to each results file
process make_matrix_and_headers {

    stageInMode 'copy'

    input:
    file(allfiles) from outfilechannel

    output:
    file("*.{txt,tsv}") into all_files_for_tar

"""
#!/usr/bin/env python3

import django
import datetime
import os
import glob

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
django.setup()
from phenotypePredictionApp.models import *

cutoff = "${params.accuracy_cutoff}"
now = datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")
HEADER = "Model_name\\tVerdict\\tProbability\\tBalanced_accuracy\\n"

modelvec = list(set([x.model_name for x in model.objects.filter()]))
modelvec = sorted(modelvec)

bin_dict = {}
countdict = {x: {"YES": 0, "NO": 0, "N/A": 0} for x in modelvec}
resultmat = []

for name in glob.glob("*.results"):

    # extract info for matrix writing
    with open(os.path.join(os.getcwd(), name), "r") as binfile:
        binname = name.replace(".results", "")
        bin_dict[binname] = {}
        for line in binfile:
            sline = line.split()
            modelname, verdict = sline[0], sline[1]
            bin_dict[binname][modelname] = verdict
            countdict[modelname][verdict] += 1
    
        # sort results file and prepend header
        with open(os.path.join(os.getcwd(), name + ".txt"), "w") as sortfile:
            binfile.seek(0, 0)
            content = []
            for line in binfile:
                content.append(line.split())
            content = sorted(content, key=lambda x: x[0])
            sortfile.write(HEADER)
            for tup in content:
                sortfile.write("\\t".join(tup))
                sortfile.write("\\n")

countlist = []
for cond in ("YES", "NO", "N/A"):
    condlist = []
    for key, vals in countdict.items():
        condlist.append(vals[cond])
    countlist.append("\\t".join([cond] + [str(x) for x in condlist]))

with open("summary_matrix.results.tsv", "w") as outfile:
    outfile.write("# phenDB\\n# Time of run: {da}\\n# Accuracy cut-off: {co}\\n".format(da=now, co=cutoff))
    outfile.write("#\\nSummary of models:\\n")
    outfile.write("\\t".join([" "] + modelvec))
    outfile.write("\\n")
    for line in countlist:
        outfile.write(line)
        outfile.write("\\n")

with open("per_bin_matrix.results.tsv", "w") as outfile2:
    outfile2.write("# phenDB\\n# Time of run: {da}\\n# Accuracy cut-off: {co}\\n".format(da=now, co=cutoff))
    outfile2.write("\\n#Results per bin:\\n")
    outfile2.write("\\t".join([" "] + modelvec))
    outfile2.write("\\n")
    for item in bin_dict.keys():
        resultlist = []
        for modelname in modelvec:
            try:
                resultlist.append(bin_dict[item][modelname])
            except KeyError:
                resultlist.append("NC")
        outfile2.write("\\t".join([item] + resultlist))
        outfile2.write("\\n")
"""
}


process zip_results {

    tag { jobname }

    stageInMode 'copy'  //actually copy in the results so we not only tar symlinks

    input:
    file(allfiles) from all_files_for_tar
    file(errorfile)

    output:
    file("${jobname}.zip") into zip_to_db

    script:
    """
    mkdir -p ${jobname}/summaries
    cp ${errorfile} ${jobname}/summaries/input_errors.log
    mv *.results.tsv ${jobname}/summaries
    mv *.results.txt ${jobname}
    zip -r ${jobname}.zip ./${jobname}
    rm -rf ${jobname}
    """
}

db_written = picaout_db_write.collectFile() { item ->
    [ "${item[1]}.results", "${item[2]}\t${item[3]}\t${item[4]}" ]  // use md5sum as filename
}


//TODO: add pica_pvalue to model_results
process write_pica_result_to_db {

    //errorStrategy 'ignore'

    input:
    file(mdsum_file) from db_written

    output:
    stdout exo

    script:
// language=Python
"""
#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import UploadedFile, bin, model, result_model

# get bin from db
try:
    parentbin = bin.objects.get(md5sum="${mdsum_file.getBaseName()}")
except ObjectDoesNotExist:
    sys.exit("Bin not found.")

conditions = []
with open("${mdsum_file}", "r") as picaresults:
    for line in picaresults:
        conditions.append(line.split())

modelresultlist=[]
for result in conditions:  # result = [modelname, verdict, pica_p_val, balanced_accuracy]
    try:
        get_bool = {"YES": True, "NO": False, "N/A": None}
        get_pica_pval = float(result[2]) if type(result[2]) != str else float(0)
        
        boolean_verdict = get_bool[result[1]]
        #get model from db
        try:
            this_model=model.objects.filter(model_name=result[0]).latest('model_train_date')
        except ObjectDoesNotExist:
            sys.exit("Current Model for this result not found.")
        
        modelresult = result_model(verdict=boolean_verdict,
                                   accuracy=float(result[3]),
                                   pica_pval=get_pica_pval,
                                   bin=parentbin,
                                   model=this_model
                                   )
        modelresultlist.append(modelresult)
    except (IntegrityError, ) as e:
        sys.exit(e)
result_model.objects.bulk_create(modelresultlist)
"""
}

process write_zip_to_db {

    input:
    file(zip) from zip_to_db

    script:
    errors_occurred = errorfile.isEmpty() ? "False" : "True"
// language=Python
"""
#!/usr/bin/env python3

import django
import sys
import os
from django.core.exceptions import ObjectDoesNotExist
from django.db import IntegrityError
from django.core.files import File
os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

django.setup()
from phenotypePredictionApp.models import UploadedFile

try:
    obj = UploadedFile.objects.filter(key='${jobname}')
    obj.update(errors = ${errors_occurred})
    file = open('${zip}', 'rb')
    djangoFile = File(file)
    obj[0].fileOutput.save('${jobname}.zip', djangoFile, save="True")
    obj.update(job_status = '100')
except IntegrityError:
    sys.exit("Exited with integrity error upon adding results to database.")
"""
}

workflow.onComplete {
    println "picaPipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}

workflow.onError {
    println "Pipeline has failed fatally with the error message: \n\n$workflow.errorMessage\n\n"
    println "Writing error report to directory ${outdir}/logs..."
    fatal_error_file = file("${outdir}/logs/errorReport.log")
    fatal_error_file.text = workflow.errorReport
    // maybe have to set pythonpath too
    "python3 scripts/set_pipeline_error.py ${jobname}".execute()
}
