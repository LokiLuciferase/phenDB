#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 6/12/18.
#

import os
import sys

from redis import Redis
from rq import Queue

from enqueue_job import update_taxonomy
from phenotypePredictionApp.variables import PHENDB_QUEUE, PHENDB_BASEDIR

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
sys.path.append(PHENDB_BASEDIR)

# enqueue a call to update_taxonomy() into redis queue
# which updates kronaTools taxonomy.tab and drops and rebuilds Taxon table in database

q = Queue(PHENDB_QUEUE, connection=Redis())
q.enqueue_call(func=update_taxonomy,
               args=(PHENDB_BASEDIR,),
               timeout='10m',
               ttl='10m',
               job_id='db_maintenance_update_taxonomy',
               at_front=True
               )
