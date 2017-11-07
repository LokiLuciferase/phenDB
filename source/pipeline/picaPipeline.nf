// #!/usr/bin/env nextflow

// define static variables
file(params.workdir).mkdir()
outdir = "$params.workdir/${file(params.inputfolder).getBaseName()}_results/"
file(outdir).mkdir()
models = file(params.modelfolder).listFiles()
input_files = Channel.fromPath("${params.inputfolder}/*.fasta")
input_gzipfiles = Channel.fromPath("${params.inputfolder}/*.tar.gz")
all_input_files = Channel.fromPath("${params.inputfolder}/*")
hmmdb = file(params.hmmdb)

log.info"""
    ##################################################################################
    
    PhenDB Pipeline started.

    Input folder containing .fasta genomes (--inputfolder): $params.inputfolder
    Output directory: $outdir
    
    Disabled compute nodes (for hmmer computation) (--omit_nodes): $params.omit_nodes
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
    stdout into fasta_errorchannel

    script:
    binname = item.getBaseName()
    """
    #!/usr/bin/env python
    from Bio import SeqIO
    
    with open("sanitychecked.fasta","w") as outfile:
        for read in SeqIO.parse("${item}", "fasta"):
            SeqIO.write(read, outfile, "fasta")        
    """

}
// Prints all errors of the fasta_sanity_check to file
fasta_errorchannel.subscribe {errorfile.append(it)}

// Print the number of fasta files to a file for progress display
all_fasta_input_files2.count().subscribe {fastafilecount.text=it} //LL: beautiful!


// call prodigal for every sample in parallel
// output each result as a set of the sample id and the path to the prodigal outfile
process prodigal {

    module "prodigal"
    memory = "2 GB"

    input:
    set val(binname), file(item) from fasta_sanitycheck_out

    output:
    set val(binname), file("prodigalout.faa") into prodigalout

    script:
    """
    prodigal -i ${item} -a prodigalout.faa > /dev/null
    """
}

// call hmmer daemon for every sample in series
// very un-nextflow, but necessary: after each completed round, the ID is written to global file progressbar_hmmer
// ((number of lines in file/number of samples)*100)) == overall progress of hmmer
process hmmer {

    maxForks 1  //do not parallelize!
    module "hmmer"

    input:
    set val(binname), file(item) from prodigalout

    output:
    set val(binname), file("hmmer.out"), file(item) into hmmerout

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

    module "muscle"

    input:
    set val(binname), file(hmmeritem), file(prodigalitem) from hmmerout

    output:
    set val(binname), file(hmmeritem), file(prodigalitem), file("complecontaitem.txt") into complecontaout

    """
    python2 $params.compleconta_path $prodigalitem $hmmeritem | tail -1 > complecontaitem.txt
    """
}

// compute accuracy from compleconta output and model intrinsics (once for each model).
process accuracy {

    memory = '10 MB'
    errorStrategy 'ignore'  //model files not yet complete, TODO: remove this!!!!

    input:
    set val(binname), file(hmmeritem), file(prodigalitem), file(complecontaitem) from complecontaout
    each model from models

    output:
    set val(binname), val(model), file(hmmeritem), file(prodigalitem), file(complecontaitem), stdout into accuracyout

    script:
    RULEBOOK = model.getBaseName()
    ACCURCYFILE = "$model/${RULEBOOK}.accuracy"
    """
    python2 $params.balanced_accuracy_path $ACCURCYFILE $complecontaitem
    """
}

// call pica for every sample for every condition in parallel
process pica {

    memory = '500 MB'
    errorStrategy 'ignore'  //model files not yet complete, TODO: remove this!!!!

    input:
    set val(binname), val(model), file(hmmeritem), file(prodigalitem), file(complecontaitem), val(accuracy) from accuracyout

    output:
    set val(binname), val(RULEBOOK), stdout, val(accuracy) into picaout  //print decision on stdout, and put stdout into return set

    script:
    RULEBOOK = model.getBaseName()
    TEST_MODEL = "$model/${RULEBOOK}.rules"
    float accuracy_cutoff = params.accuracy_cutoff as float
    float accuracy_float = accuracy as float

    if (accuracy_float >= accuracy_cutoff) {
    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f1 $hmmeritem | tr "\\n" "\\t" >> tempfile.tmp
    python2 $params.pica_path -m $TEST_MODEL -t $RULEBOOK -s tempfile.tmp > picaout.result
    echo -n \$(cat picaout.result | cut -f2 | tail -n1)
    """
    }

    else {
    """
    echo -ne "${binname}\t" > tempfile.tmp
    cut -f1 $hmmeritem | tr "\\n" "\\t" >> tempfile.tmp
    echo "N/A" > picaout.result
    echo -n \$(cat picaout.result | cut -f2 | tail -n1)
    """
    }
}

// merge all results into a file called $id.results and move each file to results folder.
picaout.collectFile() { item ->
    [ "${item[0]}.results", "${item[1]} ${item[2]} ${item[3]}" ]
}
.subscribe { it.copyTo(outdir) }



workflow.onComplete {
    println "picaPipeline completed at: $workflow.complete"
    println "Execution status: ${ workflow.success ? 'OK' : 'failed' }"
}
