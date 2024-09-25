#!/bin/bash


#SBATCH --job-name=SMK          # Job name

#SBATCH --output=FVS%j.log # Standard output and error log

#SBATCH --nodes=1                  # Request one nodes

#SBATCH --ntasks-per-node=1        # Request 1 tasks per node

#SBATCH --cpus-per-task=16

#SBATCH --time=72:00:00            # Time limit hrs:min:sec (e.g., 1 hour)

#SBATCH --mem=48GB # RAM per node #

#SBATCH --partition=cpu(all)

# Load any necessary modules or environment settings here

# Your job commands go below this line

#ulimit -n 2048
snakemake -s findviralstrains.smk --configfile config_files/no_ref_test.yml --cores 1 --rerun-incomplete
