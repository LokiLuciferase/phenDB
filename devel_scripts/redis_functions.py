#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 16/01/2018.
#
from redis import Redis
from rq import Queue, get_current_job

from time import sleep
from pprint import pprint

# this is only here because python-rq cannot import the actual worker job from same .py file.
def mockjob():
    sleep(3)
    return True

def failedjob():
    sleep(3)
    raise RuntimeError("pi equals three")

# TODO: these functions could be used in the django server to query the redis queue
# returns integer of jobs in queue
def query_queue_length():
    que = Queue('phenDB', connection=Redis())
    return len(que)

# returns job objects of queue themselves
def query_queue_jobs():
    que = Queue('phenDB', connection=Redis())
    return que.jobs

# returns a tuple: (str(job_id), datetime(queue_time), tuple(arguments), dict(meta))
def query_queue_meta():
    que = Queue('phenDB', connection=Redis())
    jobs = que.jobs
    metas = [(x.id, x.enqueued_at, x.args, x.meta) for x in jobs]
    return metas

# print jobs from failed queue with meta data
def examine_failed_jobs():
    fque = Queue('failed', connection=Redis())
    jobs = fque.jobs
    metas = [(x.id, x.enqueued_at, x.args, x.meta) for x in jobs]
    return metas

# delete jobs from failed queue
def drop_failed_jobs():
    fque = Queue('failed', connection=Redis())
    fque.empty()
    return True

# returns a tuple: (int(position_in_queue), int(len(queue))) or a tuple (-1, int(len(queue))) if key was not found
def get_current_position(keyname):
    que = Queue('phenDB', connection=Redis())
    jobs = list(que.jobs)
    for index, job in enumerate(jobs):
        if job.id == keyname:
            return (index, len(jobs))
    return (-1, len(jobs))


if __name__ == "__main__":

    print("Queued jobs: {ql}".format(ql=query_queue_length()))
    print("\nQueued jobs: ")
    pprint(query_queue_jobs())
    print("\nQueued jobs plus meta info: ")
    pprint(query_queue_meta())
    print("\nFailed jobs plus meta info: ")
    pprint(examine_failed_jobs())

    #print("\nCancelling failed jobs...")
    #drop_failed_jobs()
    #pprint(examine_failed_jobs())
    print("\nGetting specific job position:")
    print(get_current_position("a_real_job_key"))
