#!/usr/bin/env python3
import subprocess
import os
import os.path
import threading

#infolder = "/home/phen_work/test_samplesets/abc_folder"
infolder = "/home/phen_work/test_samplesets/test/C2406.fasta"
#infolder = "/home/phen_work/test_samplesets/dereplicated_genomes"

runscript_path = "/apps/phenDB_devel_PP/phenDB/source/pipeline/run_picaPipeline.sh"
pipeline_path = "/apps/phenDB_devel_PP/phenDB/source/pipeline/picaPipeline.nf"
above_workfolder = "/apps/phenDB_devel_PP/pipelineresults/resultFiles"

pica_cutoff = "0.5"
node_offs = ""

os.environ["DJANGO_SETTINGS_MODULE"] = "phenotypePrediction.settings"
os.environ["PYTHONPATH"] = "/apps/phenDB_devel_PP/phenDB/source/web_server:$PYTHONPATH"

# create workfolder
#outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn="dereplicated_genomes"))
#outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn="abc2"))
outfolder = os.path.join(above_workfolder, "{jn}_results".format(jn="C2406"))

os.makedirs(outfolder, exist_ok=True)

# create log folder
logfolder = os.path.join(outfolder, "logs")
os.makedirs(logfolder)

# call minimal runscript which performs nohup and output rerouting
subprocess.run([runscript_path,
                pipeline_path,
                infolder,
                outfolder,
                pica_cutoff,
                node_offs], check=True)