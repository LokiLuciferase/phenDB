#
# Created by Lukas Lüftinger on 6/15/18.
#
from phenotypePredictionApp.variables import PHENDB_QUEUE, PHENDB_BASEDIR
from enqueue_job import update_taxonomy

update_taxonomy(PHENDB_BASEDIR)
