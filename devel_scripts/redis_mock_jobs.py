#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 16/01/2018.
#

from redis import Redis
from rq import Queue
from redis_functions import mockjob, failedjob
from uuid import uuid4

# enqueues fake jobs
def queue_mock_jobs(number):

    que = Queue('phenDB', connection=Redis())
    for i in range(number):
        job = que.enqueue_call(func=mockjob,
                               job_id=str(uuid4()),
                               timeout='5h',
                               ttl='5h')
        # any desired property of the job can be set by adding it to the meta dict
        job.meta['ip'] = "192.168.8.1"
        job.meta['mail'] = "test@gmail.com"
        job.meta['filesize'] = 300
        job.save()


# enqueues fake jobs which will fail
def queue_failed_jobs(number):

    que = Queue('phenDB', connection=Redis())
    for i in range(number):
        job = que.enqueue_call(func=failedjob,
                               job_id=str(uuid4()),
                               timeout='5h',
                               ttl='5h')
        job.meta['ip'] = "10.0.0.1"
        job.meta['mail'] = "test@gmx.com"
        job.meta['filesize'] = 99999
        job.save()

queue_mock_jobs(10)
queue_failed_jobs(10)
queue_mock_jobs(10)