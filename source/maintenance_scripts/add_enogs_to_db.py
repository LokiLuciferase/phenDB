import click
import django
import gzip

django.setup()
from phenotypePredictionApp.models import Enog


@click.command()
@click.argument("f", type=click.Path(dir_okay=False, exists=True))
@click.option("--drop", is_flag=True)
def add_enogs_to_db(f, drop):
    if drop:
        Enog.objects.all().delete()
    enogs = []
    with gzip.open(f, mode="rt") as fin:
        for i, line in enumerate(fin.readlines()):
            line = line.split("\t")
            enogs.append(Enog(enog_name=line[0], enog_descr=line[1].rstrip()))
            print(f"Adding ENOG #{i}", end="\r")
    print("\nWriting ENOGs to DB...")
    Enog.objects.bulk_create(enogs)


if __name__ == "__main__":
    add_enogs_to_db()
