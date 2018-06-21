#
# Created by Lukas Lüftinger on 6/15/18.
#
from phenotypePredictionApp.variables import PHENDB_QUEUE, PHENDB_BASEDIR
from businessLogic.enqueue_job import update_taxonomy

# call the function update_taxonomy from the queued scripts section
update_taxonomy(PHENDB_BASEDIR)
