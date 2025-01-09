#!/bin/bash


#SBATCH --job-name=my_job          # Job name

#SBATCH --output=job_output_%j.log # Standard output and error log

#SBATCH --partition=cpu(all)       # Partition (queue) name

#SBATCH --nodes=2                  # Request two nodes

#SBATCH --ntasks-per-node=3        # Request 3 tasks per node

#SBATCH --time=01:00:00            # Time limit hrs:min:sec (e.g., 1 hour)


# Load any necessary modules or environment settings here


# Your job commands go below this line

snakemake -s findviralstrains.smk --configfile config_files/no_ref_test.yml --cores 2 --rerun-incomplete
