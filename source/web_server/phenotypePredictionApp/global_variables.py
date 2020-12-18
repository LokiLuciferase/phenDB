import os

DEFAULT_VALUES = {
    "balanced_accuracy_cutoff": 0.75,
    "prediction_confidence_cutoff": 0.6,
    "example_file_path": f"{os.environ['PHENDB_DATA_DIR']}/examples/Escherichia_coli_K12_MG1655.fna.gz",
    "example_file_alert": "Escherichia coli str. K-12 substr. MG1655 (Refseq accession GCF 000005845.2) will be used in this example",
}
