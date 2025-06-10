import os


config["run_location"] = os.getcwd() if config["run_location"] == "." else config["run_location"]

# DO NOT PUT SPACES IN THIS FILE! It breaks snakemake and gives awful errors - McKayl #
#################
##   GLOBALS   ##
#################
# Main config settings
ANALYSIS = config["analysis_ID"]
READ_DIR = config["read_dir"]
OUTPUT_DIR = config["output_dir"]
SEQUENCER = config["sequencer"]
READ_PURGE_PERCENT = config["read_purge_percent"]
DECOMP_TIME_LIMIT = config["decomp_time_limit"]
GUROBI_THREADS = config["gurobi_threads"]
RUN_LOCATION = os.getcwd() if config["run_location"] == "." else config["run_location"]
PRUNE_COUNT = config["prune"]
VISUALIZE = config["visualize"]
###############
##   SETUP   ##
###############


module pipeline1:
    snakefile: "findviralstrains.smk"
    config: config
 

module pipeline2:
    snakefile: "findviralstrains_2.smk"
    config: config

# Use rules from both modules with proper prefixes
use rule * from pipeline1 
use rule * from pipeline2 
