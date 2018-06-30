#!/usr/bin/env python3
import os
import sys
import threading

from redis import Redis
from rq import Queue


class StartProcessThread(threading.Thread):
    def __init__(self, keyname, req_balac):
        threading.Thread.__init__(self)
        self.keyname = keyname
        self.req_balac = req_balac

    def run(self):
        from phenotypePredictionApp.variables import PHENDB_BASEDIR, PHENDB_QUEUE, PHENDB_DEBUG
        from businessLogic.enqueue_job import phenDB_enqueue

        ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
        os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
        sys.path.append(ppath)

        ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
        infolder_base = os.path.join(PHENDB_BASEDIR, "data/uploads")
        pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/picaPipeline.nf")
        above_workfolder = os.path.join(PHENDB_BASEDIR, "data/results")
        infolder = os.path.join(infolder_base, self.keyname)

        pica_cutoff = float(self.req_balac)
        node_offs = ""

        # create workfolder and logfolder
        outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn=self.keyname))
        logfolder = os.path.join(outfolder, "logs")
        os.makedirs(outfolder, exist_ok=True)
        os.makedirs(logfolder)

        # add the function call to the redis queue
        q = Queue(PHENDB_QUEUE, connection=Redis())
        q.enqueue_call(func=phenDB_enqueue,
                       args=(ppath, pipeline_path, infolder, outfolder, pica_cutoff, node_offs),
                       timeout='48h',
                       ttl='48h',
                       job_id=self.keyname
                       )
        # return pipeline_job
        # here we could return len(q). or fetch it somewhere else. We could also set errors in the DB.
