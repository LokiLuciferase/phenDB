#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 08/05/2018.
#
import os

from phenotypePredictionApp.variables import *
from enqueue_job import phenDB_enqueue
from redis import Redis
from rq import Queue


ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
infolder = os.path.join(PHENDB_BASEDIR, "data/uploads/PHENDB_PRECALC")
outfolder = os.path.join(PHENDB_BASEDIR, "data/results/PHENDB_PRECALC_results")
pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/picaPipeline.nf")

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
os.environ["PYTHONPATH"] = ppath

# create workfolder
os.makedirs(outfolder, exist_ok=True)

# create log folder
logfolder = os.path.join(outfolder, "logs")
os.makedirs(logfolder, exist_ok=True)

print("Submitting precalculation job. Bins in folder {inf} will be added to the database.".format(inf=infolder))

q = Queue(PHENDB_QUEUE, connection=Redis())
pipeline_errorcode = q.enqueue_call(func=phenDB_enqueue,
                                    args=(ppath, pipeline_path, infolder, outfolder, 0.5, ""),
                                    timeout='72h',
                                    ttl='72h',
                                    job_id="PHENDB_PRECALC"
                                    )
