import os


# Set run location
config["run_location"] = os.getcwd() if config["run_location"] == "." else config["run_location"]


# Include pipelines with modern module syntax
module pipeline1:
    snakefile: "findviralstrains.smk"
    config: config

module pipeline2:
    snakefile: "findviralstrains_2.smk"
    config: config

# Use rules from both modules
use rule * from pipeline1 as other_*
use rule * from pipeline2 as other_*

