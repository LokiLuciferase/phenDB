import sys
import os
import os.path
import subprocess
from pathlib import Path

import click
import django


django.setup()
from phenotypePredictionApp.models import Taxon


@click.command()
@click.argument(
    'taxonomy_dir',
    default="/apps/miniconda3/opt/krona/taxonomy",
    type=click.Path(file_okay=False)
)
@click.option('--drop', is_flag=True)
def update_taxonomy(taxonomy_dir, drop):
    taxonomy_dir = Path(str(taxonomy_dir))
    taxonomy_file = taxonomy_dir/'taxonomy.tab'
    print("Updating KronaTools Taxonomy DB using ktUpdateTaxonomy.sh...")
    update_call = subprocess.run(
        ["ktUpdateTaxonomy.sh", taxonomy_dir],
        check=True
    )
    if not update_call.returncode == 0:
        raise RuntimeError("Updating of taxonomy.tab exited with error code. Aborting.")
    print("Done.")

    if not os.path.exists(taxonomy_file):
        raise RuntimeError("taxonomy.tab was not found. Cannot update Taxonomy.")

    print("Reading rows from updated taxonomy file...")
    taxonomy_entries = []
    counter = 0
    with open(taxonomy_file, "r") as taxfile:
        for line in taxfile:
            sys.stdout.write("Reading line {i}               \r".format(i=counter))
            sys.stdout.flush()
            counter += 1
            tax_id, other, parent, rank, namestring = line.strip().split("\t")
            new_taxon_tup = (tax_id, rank, namestring)
            taxonomy_entries.append(new_taxon_tup)
    print("\nDone.")

    if len(taxonomy_entries) < 1700000:
        raise RuntimeError("Something went wrong during database read-in. Aborting.")

    # drop all rows from old taxonomy DB
    if drop:
        print('Truncating taxonomy DB.')
        Taxon.objects.all().delete()
    print("Updating taxonomy DB...")
    taxa = [Taxon(tax_id=e[0], taxon_rank=e[1], taxon_name=e[2]) for e in taxonomy_entries]
    Taxon.objects.bulk_create(taxa)
    print("Finished updating Taxonomy Table.")


if __name__ == '__main__':
    update_taxonomy()
