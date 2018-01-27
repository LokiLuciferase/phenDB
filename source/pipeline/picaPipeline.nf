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

all_input_files = Channel.fromPath("${params.inputfolder}/*").filter{!it.name.endsWith('.tar.gz') && !it.name.endsWith('.zip')}

//todo: job status update is incorrect if some bins are already known

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

// initialize error file for "sanity" errors (eg. corrupt fasta files)
errorfile= file("$outdir/logs/sanity_errors.log")
errorfile.text=""

// unzip tar.gz files
process tgz {

    scratch true

    input:
    file(tarfile) from input_gzipfiles

    output:
    //file("${tarfile.getSimpleName()}/*.fasta") into tgz_unraveled_fasta
    file("${tarfile.getSimpleName()}/*") into tgz_unraveled_all

    script:
    outfolder = tarfile.getSimpleName()
    """
    mkdir ${outfolder}
    mkdir tempfolder
    tar -xf ${tarfile} --directory tempfolder
    cd tempfolder
    mv \$(find . -type f) ../${outfolder}/.
    """
}

// unzip .zip files
process unzip {

    scratch true

    input:
    file(zipfile) from input_barezipfiles

    output:
    //file("${zipfile.getSimpleName()}/*/*.fasta") into zip_unraveled_fasta
    file("${zipfile.getSimpleName()}/*") into zip_unraveled_all

    script:
    outfolder = zipfile.getSimpleName()
    """
    mkdir ${outfolder}
    mkdir tempfolder
    unzip ${zipfile} -d tempfolder
    cd tempfolder
    mv \$(find . -type f) ../${outfolder}/.
    """
}


// combine raw fasta files and those extracted from archive files
//all_fasta_input_files = input_files.mix(tgz_unraveled_fasta.flatten(), zip_unraveled_fasta.flatten())
truly_all_input_files = all_input_files.mix(tgz_unraveled_all.flatten(), zip_unraveled_all.flatten())

// Passes every fasta file through Biopythons SeqIO to check for corrupted files
process fasta_sanity_check {
//todo: output also sequences-only for md5sum?
    tag { binname }
    errorStrategy 'ignore'
    maxForks 5
    scratch true

    input:
    file(item) from truly_all_input_files

    output:
    set val(binname), file("sanitychecked.fasta") into sanity_check_for_continue, sanity_check_for_count

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
            with open("${errorfile}", "a") as myfile:
                myfile.write("WARNING: ${binname} is a compressed directory, not a file. It is dropped from further analysis; Nested directories cannot be analyzed! \\n\\n")
            sys.exit("${binname} is a compressed directory, not a file. Aborting.")
             

elif ("${item}".endswith("tar")):
    os.makedirs("extr")
    tar = tarfile.open("${item}", "r:")
    for tarinfo in tar:
        if tarinfo.isreg():
            tar.extractall(path="extr")
            tar.close()
            inputfasta="extr/"+str(os.listdir("extr")[0])
    
        else:
            with open("${errorfile}", "a") as myfile:
                myfile.write("WARNING: ${binname} is a compressed directory, not a file. It is dropped from further analysis; Nested directories cannot be analyzed! \\n\\n")
            sys.exit("${binname} is a compressed directory, not a file. Aborting.")
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

nr_of_files = sanity_check_for_count.count()

// sanity_check_for_count is later used in update_job_completeness to count the total nr. of bins
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
// todo: job completeness update
process add_bin_to_db {

    input:
    set val(binname), val(mdsum), file(item) from md5_out
    val(nr_of_files)

    output:
    set val(binname), val(mdsum), file(item), stdout, file("reconstructed_hmmer_file.txt"), file("reconstructed_compleconta_file.txt") into bin_is_in_db, bin_is_not_in_db

    tag { binname }

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
    parentjob.total_bins = ${nr_of_files}
    parentjob.save()

except ObjectDoesNotExist:
    sys.exit("Job not found.")

# if the bin exists, calculate it's hmmer and compleconta file from the db entries
try:
    thisbin=bin.objects.get(md5sum="${mdsum}")
    with open("reconstructed_hmmer_file.txt", "w") as hmmer:
        for entry in result_enog.objects.filter(bin=thisbin):
            hmmer.write("dummy\\t"+entry.enog.enog_name+"\\tdummy\\n")
        
    with open("reconstructed_compleconta_file.txt", "w") as complecon:
        complecon.write(str(thisbin.comple)+"\\t"+str(thisbin.conta))

    print("NO", end='')
    
    # TODO: also update the job completeness here, because then hmmer wont be called
    
except ObjectDoesNotExist:
    # write impossible Nr. for comple and conta that would cause an error during get_accuracy if not overwritten:
    thisbin = bin(bin_name="${binname}", md5sum="${mdsum}", comple=2, conta=2)
    thisbin.save()
    print("YES", end='')
    
    #if it not in the db yet, just create dummy files
    with open("reconstructed_hmmer_file.txt", "w") as hmmer:
        hmmer.write("dummy")
    with open("reconstructed_compleconta_file.txt", "w") as complecon:
        complecon.write("dummy")



try:
    assoc= bins_in_UploadedFile(bin=thisbin, UploadedFile=parentjob)
    assoc.save()    
except IntegrityError:
    sys.exit("Cannot add bin to db: An identical file (same md5sum) from the same job is already in the db. Please remove duplicate files from your input!")
"""

}

// call prodigal for every sample in parallel
// output each result as a set of the sample id and the path to the prodigal outfile
// only execute prodigal if the been has not already been in our db
process prodigal {

    tag { binname }
    maxForks 5

    memory = "450 MB"

    input:
    set val(binname), val(mdsum), file(item), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_hmmer) from bin_is_not_in_db

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
    set val(binname), val(mdsum), file(item), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_hmmer) from bin_is_in_db
    each model from models

    output:
    set val(binname), val(mdsum), val(model), stdout ,file(reconstr_hmmer), file(reconst_hmmer) into determined_if_model_recalculation_needed

    when:
    calc_bin_or_not == "NO" && model.isDirectory()

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
process uptodate_model_to_targz {

    input:
    set val(binname), val(mdsum), val(model), val(calc_model_or_not), file(reconstr_hmmer), file(reconst_hmmer) from dont_calc_model

    output:
    set val(binname), val(mdsum), val(RULEBOOK), stdout into new_model_to_targz_out
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

    picaresult=result_model.objects.get(model=model.objects.filter(model_name="${RULEBOOK}").latest('model_train_date'), bin=bin.objects.get(md5sum="${mdsum}"))     
    
    try:
        if picaresult.verdict == True:
            verdict="YES"
        elif picaresult.verdict == None:
            verdict = "N/A"
        else:
            verdict= "NO" 
    except:
            verdict="N/A"
    pica_pval= "N/A" if picaresult.pica_pval == 0.0 else str(picaresult.pica_pval)  
    accuracy=str(picaresult.accuracy)      
    write_this=verdict+" "+pica_pval+"\\t"+accuracy    
    print(write_this, end="")        
    
    """
}

// split output of last process on tab to get seperate verdict and accuracy fields
picaout_from_new_model = new_model_to_targz_out.map { l ->

    splitted=l[4].split("\t")
    verdict=splitted[0]
    accuracy=splitted[1] as float
    //accuracy+="\n"
    return [ l[0], l[1], l[2], verdict, accuracy ]

}

// if the results in our db for this bin and model are outdated, pica needs to be called
// this just reformats the input and just calls pica if the value is YES

accuracy_in_from_old_model = calc_model.filter{ it[3] == "YES" }.map{ l -> [l[0], l[1], l[2], l[4], l[5]] }


// call hmmer daemon for all bins for which prodigal was run (i.e. bins that have not yet been in the db)
process hmmer {

    tag { binname }
    scratch true

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
    val(nr_of_files)

    output:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) into job_updated_out

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

try:
    parentjob = UploadedFile.objects.get(key="${jobname}")
    current_finished = parentjob.finished_bins
    total_jobs = parentjob.total_bins
    if current_finished < total_jobs - 1:
        parentjob.finished_bins = int(current_finished) + 1
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

// determine if the uploaded bin is a bacterium. add verdict to set
// later on, use verdict to omit certain models for archaea
// TODO: implement this, right now all are treated as bacteria
process pica_bacteria {

    tag { binname }
    scratch true

    memory = '500 MB'

    input:
    set val(binname), val(mdsum), file(hmmeritem), val(complecontaitem) from complecontaout

    output:
    set val(binname), val(mdsum), file(hmmeritem), val(complecontaitem), stdout into complecontaout_for_call_accuracy, complecontaout_for_hmmerresults_to_db  // "YES"= is bacterium

    script:
    RULEBOOK = file(params.archaea_model_path).getBaseName()
    TEST_MODEL = "${params.archaea_model_path}/${RULEBOOK}.rules"
    """
    echo -ne "${binname}\\t" > tempfile.tmp
    cut -f2 $hmmeritem | tr "\\n" "\\t" >> tempfile.tmp
    test.py -m $TEST_MODEL -t $RULEBOOK -s tempfile.tmp > picaout.result
    is_archaea=\$(cat picaout.result | tail -n1 | cut -f2)
    if [[ \$is_archaea = "YES" ]]; then
        echo -n "NO"
    else
        echo -n "YES"
    fi
    """
}


process write_hmmer_results_to_db {

    tag { binname }

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(complecontaitem), val(is_bacterium) from complecontaout_for_hmmerresults_to_db

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
            #todo: should that not throw an error?
        except IntegrityError:
            print("Skipping due to IntegrityError.")
            continue     
            

    try:
        result_enog.objects.bulk_create(enog_res_objectlist)        
    except IntegrityError:
        print("Could not add enogs of bin ${binname} to the db.")  
    
"""
}

// instantiates an item for each compleconta output item and each model,
// drops all channel items which are not either bacteria
// or archea but the assorted model is not contained in the blacklist omit_in_archaea
accuracy_in = complecontaout_for_call_accuracy
                .combine(models)
                .filter{it[4] == "YES" || !(params.omit_in_archaea.contains(it[5].getBaseName()))}
                .map { l -> [ l[0], l[1], l[5], l[2], l[3] ]}

// compute accuracy from compleconta output and model intrinsics.
process accuracy {

    tag { "${binname}_${model.getBaseName()}" }

    memory = '300 MB'

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
        parentbin.comple = float(cc[0])
        parentbin.conta= float(cc[1])
        parentbin.save()


        
    except ObjectDoesNotExist:
        sys.exit("Bin not found.")
    
    #check the balanced accuracy. other metrices would be very similar to evaluate, just change the part after the last dot
    # the print statments includes a newline at the end, this is important for processing further downstream
    
    print(model_accuracies.objects.get(model=model.objects.filter(model_name="${model.getBaseName()}").latest('model_train_date'),
    comple=round_nearest(float(cc[0]),0.05), 
    conta=round_nearest(float(cc[1]),0.05)).mean_balanced_accuracy)

    """
}

// call pica for every sample for every condition in parallel
process pica {

    tag { "${binname}_${model.getBaseName()}" }
    scratch true
    memory = '800 MB'

    input:
    set val(binname), val(mdsum), val(model), file(hmmeritem), val(accuracy) from accuracyout

    output:
    set val(binname), val(mdsum), val(RULEBOOK), stdout, val(accuracy) into picaout_db_write, picaout_for_download //print decision on stdout, and put stdout into return set

    script:
    RULEBOOK = model.getBaseName()
    TEST_MODEL = "$model/${RULEBOOK}.rules"

    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f2 $hmmeritem | tr "\n" "\t" >> tempfile.tmp
    test.py -m $TEST_MODEL -t $RULEBOOK -s tempfile.tmp > picaout.result
    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3)
    """
}

// replace the verdict of pica with "N/A N/A" if the balanced accuracy is below the treshold
NA_replaced_for_download = picaout_for_download.mix(picaout_from_new_model).map { l ->

    float accuracy_cutoff = params.accuracy_cutoff as float
    float accuracy_float = l[4] as float

    if (accuracy_float >= accuracy_cutoff) {
        output = l[3]
    } else {
        output = "N/A N/A"
    }
    return [l[0], l[1], l[2], output, l[4]]
}

// merge all results into a file called $id.results.
outfilechannel = NA_replaced_for_download.collectFile() { item ->
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
// language=python
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
    scratch true
    stageInMode 'copy'  //actually copy in the results so we not only tar symlinks

    input:
    file(allfiles) from all_files_for_tar
    file(errorfile)

    output:
    file("${jobname}.zip") into zip_to_db

    script:
    """
    mkdir -p ${jobname}/summaries
    cp ${params.description_file} ${jobname}/summaries/PICA_trait_descriptions.txt
    cp ${errorfile} ${jobname}/summaries/invalid_input_files.log.txt
    mv *.results.tsv ${jobname}/summaries
    mv *.results.txt ${jobname}
    zip -r ${jobname}.zip ./${jobname}
    rm -rf ${jobname}
    """
}

db_written = picaout_db_write.collectFile() { item ->
    [ "${item[1]}.results", "${item[2]}\t${item[3]}\t${item[4]}" ]  // use md5sum as filename
}

process write_pica_result_to_db {

    scratch true

    input:
    file(mdsum_file) from db_written

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
        if result[2]=="NA" or result[2]=="N/A":
            get_pica_pval=float(0)
        else:
            get_pica_pval=float(result[2])

        
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

    scratch true

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
    totalbins = obj[0].total_bins
    obj.update(finished_bins = totalbins)
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
}
