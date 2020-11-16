# Global variables for PhenDB to easily switch between development and prod instance
import os

PHENDB_BASEDIR = os.environ['PHENDB_BASEDIR']
PHENDB_DATA_DIR = os.environ['PHENDB_DATA_DIR']
PHENDB_QUEUE = os.environ['PHENDB_QUEUE']
PHENDB_DB = os.environ['PHENDB_DB_NAME']
PHENDB_DEBUG = bool(os.environ['PHENDB_DEBUG'])
