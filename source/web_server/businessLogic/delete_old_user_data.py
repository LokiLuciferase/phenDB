#!/usr/bin/env python3
#
# Created by Lukas Lüftinger on 19/02/2018.
#
from enqueue_job import delete_user_data
from redis import Redis
from rq import Queue

# enqueue a call to delete_user_data() into the redis queue
# which deletes all user data older than days_back
days_back = 30

q = Queue('phenDB', connection=Redis())
q.enqueue_call(func=delete_user_data,
               args=(days_back,),
               timeout='10m',
               ttl='10m',
               job_id="db_maintenance_clean_old",
               at_front=True
               )
