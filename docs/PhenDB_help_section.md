### PhenDB is an automated pipeline for the prediction of microbial phenotypes based on comparative genomics.



Concept
=======


The PhenDB pipeline takes as input genomic bins from a Metagenomics experiment and makes predictions regarding the traits associated with OTUs belonging to these bins.
1. bins provided as nucleotide fasta files are used to perform gene prediction with Prodigal (Hyatt et al., 2010).
2. Hmmer (Eddy et al., 2001) then searches orthologous groups for the resulting protein sequences using the EggNOG DB (Huerta-Cepas et al., 2015).
3. Compleconta (Hyden et al., unpublished) estimates completeness, contamination and taxonomic identity of the uploaded bins given a set of marker genes,
4. Finally the machine learning framework PICA (Feldbauer et al., 2015) uses support vector machine (SVM)-based models (calculated from genomes with known phenotypes) and the list of orthologous groups of proteins present in a bin/genome to predict whether the organism possesses a particular trait, and the confidence in that prediction.  
Currently, predictions for 45 different traits are calculated.  

Please note that PhenDB is still in an early stage of development and for several models we are still working on improving the training data.
Thus, please take note of the "balanced_accuracy" values ascribed to predictions.  
Additionally, please note that PhenDB currently does not provide trait predictions for archaeal bins.


How to Use PhenDB
=================

### Input File Formats
The PhenDB pipeline takes as input genomic bins from a Metagenomics experiment, which may be provided as:

* Nucleotide fasta file (raw or gzip-compressed)
* Protein fasta file (raw or gzip-compressed)
* tar.gz or .zip archive of the above

Please note that a flat file structure is required within archives.  
The current maximum file size for upload is 1 GB, the maximum file size per bin is 30 MB.  
Duplicate sequence files (determined by file content) will be silently dropped from the analysis.  
Empty sequence files will be silently dropped from the analysis.  

### Output Files
PhenDB provides a downloadable archive named after your Job ID (i.e. the key after ../results/ in the URL). This archive contains:  
*   the folder individual_results, containing:
    *   a {bin name}.results.csv file for every valid uploaded bin/genome.
    This file contains the model names, predictions (YES/NO/NA) along with PICA probability and balanced accuracy values.
*   the folder "summaries":
    *   "summary_matrix.csv": A summary file that shows for each model how many bins/genomes were predicted as "YES", "NO" or "N/A"
    *   "per\_bin\_matrix.csv": A summary file that shows the verdict for each bin and each model as a matrix
    *   "invalid\_input\_files.log.txt": If one or more of your uploaded files were invalid (e.g. not in FASTA format), a warning will appear to check this file. If all files were correct, this file is empty.
    *   "PICA\_trait\_descriptions.txt": Contains the model names and the traits they are testing for.


### Statistical Measures
PhenDB provides two separate confidence measures associated with trait predictions.
1. PICA probability is the internal probability of class membership within the linear SVM model used by PICA. It can thus be seen as the confidence the model itself shows in its prediction.
2. Balanced Accuracy is a confidence measure computed from completeness/contamination of the uploaded bin and the model's known predictive power. It can thus be interpreted as our confidence in the predictions of the model, given the bin's completeness and contamination. Values of Balanced Accuracy may range from 0.5 (predictions are expected to be random) to 1 (high expected correlation of predicted class and actual class).

A balanced accuracy cutoff can be set when submitting to PhenDB: predictions with a balanced accuracy below the chosen cutoff value (range: 0.5 - 1) are then masked in the final result (NA instead of YES/NO).

### Additional Information
After submission, your job is queued and waits for completion of any previously submitted jobs.  
When computation starts, expect about 1-1.5 min of calculation time per bin.  
To receive a notification upon job completion, enter an email address during your submission.  
Alternatively, you may save the URL to your submission to retrieve your results later.  
Results are stored by PhenDB for 30 days - after which all user-uploaded data and associated results are deleted.


References
==========


1. Hyatt, Doug, et al. "Prodigal: prokaryotic gene recognition and translation initiation site identification." BMC bioinformatics 11.1 (2010): 119.  
2. Eddy, Sean R. "HMMER: Profile hidden Markov models for biological sequence analysis." (2001).  
3. Huerta-Cepas, Jaime, et al. "eggNOG 4.5: a hierarchical orthology framework with improved functional annotations for eukaryotic, prokaryotic and viral sequences." Nucleic acids research 44.D1 (2015): D286-D293.  
4. Hyden, Patrick, et al., unpublished (https://github.com/phyden/compleconta)
5. Feldbauer, Roman, et al. "Prediction of microbial phenotypes based on comparative genomics." BMC bioinformatics 16.14 (2015): S1.
