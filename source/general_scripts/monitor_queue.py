#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 14/01/2018.
#
from redis import Redis
from rq import Queue, get_failed_queue
from rq.worker import Worker
from pprint import pprint

PHENDB_QUEUE = "phenDB_devel_LL"
#PHENDB_QUEUE = "phenDB"

rconn = Redis()
workers = Worker.all(connection=rconn)
phendb_worker = workers[0]
q = Queue(PHENDB_QUEUE, connection=rconn)

print("Successfully executed: {se}\nFailed: {fe}\nTotal time worked: {tt}\n".format(se=phendb_worker.successful_job_count,
                                                                                    fe=phendb_worker.failed_job_count,
                                                                                    tt=phendb_worker.total_working_time))
print("Currently active workers: ")
pprint([x.name for x in workers])

print("Currently enqueued jobs: ")
pprint(q.job_ids)

print("Failed jobs:")
fque = Queue('failed', connection=rconn)
pprint(fque.job_ids)
