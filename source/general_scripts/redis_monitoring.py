#!/usr/bin/env python3

#
# Created by Lukas Lüftinger on 14/01/2018.
#

from redis import Redis
from rq import Queue, Worker, get_failed_queue
from pprint import pprint


workers = Worker.all(connection=Redis())
phendb_worker = Worker.find_by_key('rq:worker:phenDB')
q = Queue('phenDB', connection=Redis())

print("Successfully executed: {se}\nFailed: {fe}\nTotal time worked: {tt}\n".format(se=phendb_worker.successful_job_count,
                                                                                    fe=phendb_worker.failed_job_count,
                                                                                    tt=phendb_worker.total_working_time))
print("Currently active workers: ")
print(workers)

print("Currently enqueued jobs: ")
pprint(q.job_ids)

print("Failed jobs:")
pprint(get_failed_queue().job_ids)
