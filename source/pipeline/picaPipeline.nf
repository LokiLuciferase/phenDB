#!/usr/bin/env nextflow
// define static variables


def flat_increment(task_attempt, increment){
    def incr = task_attempt == 1 ? 1 : task_attempt * increment
    return incr
}

// Function to ensure that resource requirements don't go beyond a maximum limit
def check_max(obj, type) {
  if(type == 'memory'){
      try {
            if((obj as nextflow.util.MemoryUnit).compareTo(params.max_memory as nextflow.util.MemoryUnit) == 1)
              return params.max_memory as nextflow.util.MemoryUnit
            else
              return obj as nextflow.util.MemoryUnit
          } catch (all) {
                println "   ### ERROR ###   Max memory '${params.max_memory}' is not valid! Using default value: $obj"
                return obj as nextflow.util.MemoryUnit
              }
    } else if(type == 'time'){
        try {
              if(obj.compareTo(params.max_time as nextflow.util.Duration) == 1)
                return params.max_time as nextflow.util.Duration
              else
                return obj
            } catch (all) {
                  println "   ### ERROR ###   Max time '${params.max_time}' is not valid! Using default value: $obj"
                  return obj
                }
    } else if(type == 'cpus'){
        try {
              return Math.min( obj, params.max_cpus as int )
            } catch (all) {
                println "   ### ERROR ###   Max cpus '${params.max_cpus}' is not valid! Using default value: $obj"
                return obj
            }
    }
}


params.recalc = false
models = file(params.modelfolder).listFiles()
hmmdb = file(params.hmmdb)

if (params.recalc) {
    jobname = "PHENDB_PRECALC"
    input_gzipfiles = Channel.empty()
    input_barezipfiles = Channel.empty()
    all_input_files = Channel.empty()
    params.batch_no = -1
    params.total_batch_no = -1
    log.info"""Recalculating batch ${params.batch_no} of ${params.total_batch_no}.""".stripIndent()
} else {
    jobname = file(params.inputfolder).getBaseName()
    input_gzipfiles = Channel.fromPath("${params.inputfolder}/*.{tar.gz,tgz}")
    input_barezipfiles = Channel.fromPath("${params.inputfolder}/*.zip")
    all_input_files = Channel.fromPath("${params.inputfolder}/*")
                      .filter{!it.name.endsWith('.tar.gz') && !it.name.endsWith('.zip') && !it.name.endsWith('.tgz')}
}

outdir = file(params.outdir)
file(outdir).mkdir()


log.info"""
    ##################################################################################

    PhenDB Pipeline started.

    Input folder containing bin files: $params.inputfolder
    Output directory: $outdir
    Job name: $jobname

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
truly_all_input_files = all_input_files.mix(tgz_unraveled_all.flatten(), zip_unraveled_all.flatten())

// check if file names contain only ascii and if file size is OK
truly_all_input_files.map {

    if (!( it.getBaseName() ==~ /^\p{ASCII}+$/ )) {
        asciimess = "WARNING: The filename of ${it.getName()} contains non-ASCII characters. " +
                    "The file was dropped from the analysis.\n\n"
        log.info(asciimess)
        errorfile.append(asciimess)
        return
    } else if ( it.size() > params.max_bin_size ) {
        sizemess = "WARNING: The size of file ${it.getName()} exceeds the maximum processible file size. " +
                   "The file was dropped from the analysis.\n\n"
        log.info(sizemess)
        errorfile.append(sizemess)
        return
    }
    return it
}.into{ groovy_checked_for_count; groovy_checked_for_continue }


// number of files for processing: before fasta checking, since waiting for this is a bottleneck
nr_of_files = groovy_checked_for_count.count()


// Pass every fasta file through Biopythons SeqIO to check for corrupted files
process fasta_sanity_check {

    tag { binname }
    errorStrategy 'ignore'  // failing files removed from pipeline
    scratch true
    maxForks 5

    input:
    file(item) from groovy_checked_for_continue

    output:
    set val(binname), file("sanitychecked.fasta"), stdout into fasta_sanity_check_out

    script:
    binname = item.getName()

    """
    fasta_sanity_check.py ${binname} ${item} ${errorfile}
    """
}

// get primary key for bins
process md5sum {

    tag { binname }

    input:
    set val(binname), file(item), val(seqtype) from fasta_sanity_check_out

    output:
    set val(binname), stdout, file(item), val(seqtype) into md5_out

    script:
    """
    echo -n \$(md5sum ${item} | cut -f1 -d" ")
    """
}

// add each bin to database
process add_bin_to_db {

    tag { binname }
    errorStrategy 'ignore'  // if duplicate files in same job, just drop them from pipeline

    input:
    set val(binname), val(mdsum), file(item), val(seqtype) from md5_out
    val(nr_of_files)

    output:
    set val(binname), val(mdsum), file(item), val(seqtype), stdout, file("reconstructed_hmmer_file.txt"), file("reconstructed_compleconta_file.txt") into bin_is_in_db, bin_is_not_in_db

    script:
    """
    add_bin_to_db.py ${jobname} ${nr_of_files} ${mdsum} ${binname}
    """
}

// call prodigal for every sample in parallel
// output each result as a set of the sample id and the path to the prodigal outfile
// only execute prodigal if the been has not already been in our db

new_dna_bins = Channel.create()
new_protein_bins = Channel.create()
bin_is_not_in_db.choice( new_dna_bins, new_protein_bins ) { l -> l[3] == "DNA" ? 0 : 1 }

process prodigal {

    tag { binname }
    memory = "450 MB"
    maxForks 5

    input:
    set val(binname), val(mdsum), file(sanitychecked), val(seqtype), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_cc) from new_dna_bins

    output:
    set val(binname), val(mdsum), file("prodigalout.faa"), val(seqtype), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_cc) into prodigalout

    when:
    calc_bin_or_not.equals("YES")

    script:
    """
    prodigal -i ${sanitychecked} -a prodigalout.faa > /dev/null
    rm -f \$(realpath ${sanitychecked})
    rm ${sanitychecked}
    """
}

// merge uploaded protein files with our calculated prodigal results
all_protein_files = prodigalout.mix( new_protein_bins.filter{ it[4] == "YES" } )

// if the bin has already been in our db:
//TODO: this is horribly inefficient, think about batch solution
process determine_models_that_need_recalculation {

    tag { "${binname}_${model.getBaseName()}" }
    errorStrategy 'ignore'

    input:
    set val(binname), val(mdsum), file(item), val(seqtype), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_cc) from bin_is_in_db
    each model from models

    output:
    set val(binname), val(mdsum), val(model), stdout, file(reconstr_hmmer), file(reconst_cc) into calc_model
    val(mdsum) into do_not_calc_model

    when:
    calc_bin_or_not == "NO" && model.isDirectory()

    script:
    // language=Python
    """
    #!/usr/bin/env python3

    import django
    import os
    os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"

    django.setup()
    from phenotypePredictionApp.models import PicaResult, PicaResultExplanation

    #if this succeeds, then there is a result for this bin and the newest model version in our db
    try:
        PicaResult.objects.get(
            model=PicaModel.objects.filter(model_name="${model.getBaseName()}").latest('model_train_date'), 
            bin=Bin.objects.get(md5sum="${mdsum}")
        )
        if "${params.get_explanations}" == "true":
            PicaResultExplanation.objects.get(
                model=PicaModel.objects.filter(model_name="${model.getBaseName()}").latest('model_train_date'),
                bin=Bin.objects.get(md5sum="${mdsum}")
            )
        print("NO", end='')
    except:
        print("YES", end='')
    """
}

// if the results in our db for this bin and model are outdated, pica needs to be called
// this just reformats the input and just calls accuracy if the value is YES
accuracy_in_from_old_model = calc_model.filter{ it[3] == "YES" }.map{ l -> [l[0], l[1], l[2], l[4], l[5]] }


// call hmmer daemon for all bins for which prodigal was run (i.e. bins that have not yet been in the db)
process hmmer {

    tag { binname }
    maxForks 1
    time 10.m

    input:
    set val(binname), val(mdsum), file(item), val(seqtype), val(calc_bin_or_not), file(reconstr_hmmer), file(reconst_cc) from all_protein_files

    output:
    set val(binname), val(mdsum), file("hmmer.out"), file(item) into hmmerout

    script:
    if (params.annotation_strategy == 'hmmer'){
    """
    DISABLED=\$(echo ${params.omit_nodes} | sed -e 's/\\([^ ]\\+\\)/-e &/g')
    if [[ -n "\$DISABLED" ]] ; then
        HMM_DAEMONCLIENTS=\$(echo cubeb{01..30} | tr " " "\\n" | grep -v \$(echo \$DISABLED) | tr "\\n" " ")
    else
        HMM_DAEMONCLIENTS=\$(echo cubeb{01..30})
    fi
    echo \$HMM_DAEMONCLIENTS

    hmmc.py -i ${item} -d $hmmdb -s \$HMM_DAEMONCLIENTS -n 100 -q 5 -m 1 -f 120 -o hmmer.out
    """
    } else {
    """
    deepnog infer -bs 1 -nw -1 -of tsv ${item} | sed 1d | awk '\$2!=""' > hmmer.out
    """
    }
}


//increase completed job count after each hmmer job completion
process update_job_completeness {

    tag { jobname }
    maxForks 1

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
    parentjob = Job.objects.get(key="${jobname}")
    current_finished = parentjob.finished_bins
    total_jobs = parentjob.total_bins
    if current_finished < total_jobs - 1:
        parentjob.finished_bins = int(current_finished) + 1
        parentjob.save()

except ObjectDoesNotExist:
    sys.exit("Job not found.")
"""
}


// compute contamination and completeness using compleconta, and taxonomy
process compleconta {

    tag { binname }

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(prodigalitem) from job_updated_out

    output:
    set val(binname), val(mdsum), file(hmmeritem), file("complecontaitem.txt") into cc_for_accuracy, cc_to_db

    script:
    """
    set -euo pipefail
    compleconta_py3.py $prodigalitem $hmmeritem | tail -1 > complecontaitem.txt
    """
}


// write output of compleconta to database for each bin
process write_hmmer_results_to_db {

    tag { binname }

    input:
    set val(binname), val(mdsum), file(hmmeritem), file(complecontaitem) from cc_to_db

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
    parentbin = Bin.objects.get(md5sum="${mdsum}")
except ObjectDoesNotExist:
    sys.exit("Bin not found.")

# write hmmer output to HmmerResult table
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
            new_enog_res = HmmerResult(bin=parentbin, enog=Enog.objects.get(enog_name=enog_elem))
            enog_res_objectlist.append(new_enog_res)
            #todo: should that not throw an error?
        except IntegrityError:
            continue
    try:
        HmmerResult.objects.bulk_create(enog_res_objectlist)
    except IntegrityError:
        print("Could not add enogs of bin ${binname} to the db.")
"""
}


// instantiates an item for each compleconta output item and each model,
accuracy_in = cc_for_accuracy
              .combine(models)
              .map { l -> [ l[0], l[1], l[4], l[2], l[3] ]}


// compute accuracy from compleconta output and model intrinsics.
process accuracy {

    tag { "${binname}_${model.getBaseName()}" }
    memory = '300 MB'

    input:
    set val(binname), val(mdsum), val(model), file(hmmeritem), file(complecontaitem) from accuracy_in.mix(accuracy_in_from_old_model)

    output:
    set val(binname), val(mdsum), val(model), file(hmmeritem), stdout into accuracyout

    script:
    modelname = model.getBaseName()
    """
    accuracy.py ${complecontaitem} ${mdsum} ${modelname}
    """
}

//make dataset to recalculate PRECALC genomes already in database; calculate batchwise
process get_recalc_hashes {

    input:
    val params.modelfolder

    output:
    file("recalc_table.csv") into recalc_table
    file("*.out") into recalc_hmmerfiles

    when:
    params.recalc

    script:
    """
    make_recalc_output_table.py ${params.modelfolder} ${params.batch_no} ${params.total_batch_no}
    """

}
recalc_table_collated = recalc_table.splitCsv(sep: "\t")
        .map{l -> [ binname: l[0],
                    mdsum: l[1],
                    model: file(l[2]),
                    hmmeritem: file(l[3]),
                    accuracy: l[4] + "\n" ]}  //same output as given by process accuracy
//recalc_table_collated = recalc_table.splitCsv(sep: "\t")
//        .map{l ->
//        if (l[0] != "dummy_bin"){
//            return [ binname: l[0], mdsum: l[1], model: file(l[2]), hmmeritem: file(l[3]), accuracy: l[4] + "\n" ]
//        } else {
//            return
//        }
//        }  //same output as given by process accuracy

pica_in = accuracyout.mix(recalc_table_collated)

// call pica for every sample for every condition
process pica {

    tag { "${binname}_${model.getBaseName()}" }
    cpus 1
    errorStrategy { task.exitStatus in [104,134,136,137,138,139,143] ? 'retry' : 'terminate' }
    memory { check_max((params.phenotrex_memory as nextflow.util.MemoryUnit) * flat_increment(task.attempt, 0.75), 'memory') }

    input:
    set val(binname), val(mdsum), val(model), file(hmmeritem), val(accuracy) from pica_in

    output:
    set val(binname), val(mdsum), val(RULEBOOK), stdout, val(accuracy) into picaout_db_files
    set val(binname), val(mdsum), val(RULEBOOK), file("*_explanations_formatted.tsv") optional true into picaout_db_files_explanations

    script:
    explain_flag = params.get_explanations ? "--out_explain_per_sample explanations.tsv --shap_n_samples ${params.shap_n_samples} --n_max_explained_features ${params.shap_n_features}" : ''
    RULEBOOK = model.getBaseName()
    TEST_MODEL = "$model/${RULEBOOK}.svm.class"

    """
    echo "${task.memory}, ${task.cpus}" > resources.log
    echo "#feature_type:eggNOG5-tax-2" > tempfile.tmp
    echo -ne "${binname}\t" >> tempfile.tmp
    cut -f2 $hmmeritem | tr "\\n" "\\t" >> tempfile.tmp
    phenotrex predict \\
        --genotype tempfile.tmp ${explain_flag} \\
        --classifier $TEST_MODEL > picaout.result
    echo -n \$(cat picaout.result | tail -n1 | cut -f2,3)
    if [[ -f "explanations.tsv" ]]; then
        cut -f1,2,3,4,5 explanations.tsv \\
            | sed 1d \\
            | awk '{print "${RULEBOOK}", \$1, \$3, \$4, \$5}' > ${mdsum}_${RULEBOOK}_explanations_formatted.tsv
    fi
    """
}

db_write_pica_results = picaout_db_files.collectFile() { item ->
    [ "${item[1]}.results", "${item[2]}\t${item[3]}\t${item[4]}" ]  // use md5sum as filename for identification in DB
}
db_write_pica_explanations = picaout_db_files_explanations.collectFile(){
    item -> ["${item[1]}.explanations", item[3].text + "\n"]
}

// write PICA result to database
process write_pica_result_to_db {

    tag "${mdsum_file.getBaseName()}"
    scratch true

    input:
    file(mdsum_file) from db_write_pica_results

    output:
    file(mdsum_file) into resfiles_after_db_write

    script:
    """
    write_pica_result_to_db.py ${mdsum_file} ${mdsum_file.getBaseName()}
    """
}

process write_pica_explanations_to_db {

    tag "${mdsum_file.getBaseName()}"
    scratch true

    input:
    file(mdsum_file) from db_write_pica_explanations

    output:
    file(mdsum_file) into explanations_after_db_write
    
    script:
    """
    write_pica_explanations_to_db.py ${mdsum_file} ${mdsum_file.getBaseName()}
    """
}

job_mdsums_expl = explanations_after_db_write.map {item -> item.getBaseName() }.unique().collect()
job_mdsums_results = resfiles_after_db_write.map { item -> item.getBaseName() }.mix(do_not_calc_model).mix(job_mdsums_expl).unique().collect()
// access database to get results for all bin md5sums received
// apply filtering by cutoffs and model constraints
// and produce individual and summary files for download
process make_downloadable_flat_files {

    input:
    val(job_mdsums_results)

    output:
    file("*.{csv,txt}") into downloadable_files_out
    file("taxonomy_krona.tsv") into tax_for_krona

    script:
    """
    write_pipeline_output.py \\
        --dep_file ${params.pica_dependencies} \\
        --job_key ${jobname} \\
        --md5sums ${job_mdsums.join(" ")}
    """
}

// make kronaplot
process krona {

 scratch true

 input:
 file(kronain) from tax_for_krona

 output:
 file('krona_taxonomy.html') into krona_out

 script:
 """
 ktImportTaxonomy $kronain -o krona_taxonomy.html
 """
}

//zip up all downloadable results
process zip_results {

    tag { jobname }
    scratch true
    stageInMode 'copy'  //actually copy in the results so we not only tar symlinks

    input:
    file(allfiles) from downloadable_files_out.collect()
    file(kronafile) from krona_out.collect()
    file(errorfile)

    output:
    file("phendb_${jobname}.zip") into zip_to_db

    script:
    """
    mkdir -p ${jobname}/summaries
    mkdir -p ${jobname}/individual_results
    cp ${errorfile} ${jobname}/summaries/invalid_input_files.log.txt
    mv *.traits.csv ${jobname}/individual_results
    mv *.csv ${jobname}/summaries
    mv *.html ${jobname}/summaries
    mv ./${jobname} phendb_${jobname}
    zip -r phendb_${jobname}.zip ./phendb_${jobname}
    rm -rf phendb_${jobname}
    """
}

//mark job as finished and save zip to upload location
process write_zip_to_db {

    scratch true

    input:
    file(zip) from zip_to_db

    script:
    errors_occurred = errorfile.isEmpty() ? "False" : "True"
    errtype = errorfile.isEmpty() ? "" : "INPUT"
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
from phenotypePredictionApp.models import Job

try:
    obj = Job.objects.filter(key='${jobname}')
    obj.update(errors=${errors_occurred}, error_type="${errtype}")
    file = open('${zip}', 'rb')
    djangoFile = File(file)
    obj[0].fileOutput.save('phendb_${jobname}.zip', djangoFile, save="True")
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
