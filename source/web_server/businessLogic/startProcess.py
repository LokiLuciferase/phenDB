#!/usr/bin/env python3
import os
import os.path
import threading

from redis import Redis
from rq import Queue

from phenotypePredictionApp.models import UploadedFile
from phenotypePredictionApp.variables import *
from businessLogic.enqueue_job import phenDB_enqueue


class StartProcessThread(threading.Thread):
    def __init__(self, keyname, req_balac):
        threading.Thread.__init__(self)
        self.keyname = keyname
        self.req_balac = req_balac

    def run(self):

        ppath = PHENDB_BASEDIR + "/source/web_server:$PYTHONPATH"
        infolder_base = os.path.join(PHENDB_BASEDIR, "/data/uploads")
        pipeline_path = os.path.join(PHENDB_BASEDIR, "source/pipeline/picaPipeline.nf")
        above_workfolder = os.path.join(PHENDB_BASEDIR, "/data/results")
        infolder = os.path.join(infolder_base, self.keyname)

        # TODO: we should make these parameters settable from the web mask
        pica_cutoff = float(self.req_balac)
        node_offs = ""

        os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
        os.environ["PYTHONPATH"] = ppath

        # create workfolder
        outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn=self.keyname))
        os.makedirs(outfolder, exist_ok=True)

        # create log folder
        logfolder = os.path.join(outfolder, "logs")
        os.makedirs(logfolder)

        # add the function call to the redis queue
        q = Queue(PHENDB_QUEUE, connection=Redis())
        pipeline_errorcode = q.enqueue_call(func=phenDB_enqueue,
                                      args=(ppath, pipeline_path, infolder, outfolder, pica_cutoff, node_offs),
                                      timeout='48h',
                                      ttl='48h',
                                      job_id=self.keyname
                                      )
        # return pipeline_job
        # here we could return len(q). or fetch it somewhere else. We could also set errors in the DB.
