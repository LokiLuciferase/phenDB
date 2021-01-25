# Global variables for PhenDB to easily switch between development and prod instance
import os

PHENDB_BASEDIR = os.environ["PHENDB_BASEDIR"]
PHENDB_DATA_DIR = os.environ["PHENDB_DATA_DIR"]
PHENDB_QUEUE = os.environ["PHENDB_QUEUE"]
PHENDB_DB = os.environ["PHENDB_DB_NAME"]
PHENDB_DEBUG = os.environ["PHENDB_DEBUG"] == 'True'
PHENDB_DEV_EMAIL = os.environ["PHENDB_DEV_EMAIL"]