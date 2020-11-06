#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 6/12/18.
#
import os
import sys

from redis import Redis
from rq import Queue

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
PHENDB_BASEDIR = os.environ['BASEDIR']
sys.path.append(f"{PHENDB_BASEDIR}/source/web_server")


from phenotypePredictionApp.variables import PHENDB_QUEUE, PHENDB_BASEDIR
from enqueue_job import update_taxonomy

# enqueue a call to update_taxonomy() into redis queue
# which updates kronaTools taxonomy.tab and drops and rebuilds Taxon table in database

q = Queue(PHENDB_QUEUE, connection=Redis())
q.enqueue_call(func=update_taxonomy,
               args=(PHENDB_BASEDIR,),
               timeout='240m',
               ttl='240m',
               job_id='db_maintenance_update_taxonomy',
               at_front=True
               )
