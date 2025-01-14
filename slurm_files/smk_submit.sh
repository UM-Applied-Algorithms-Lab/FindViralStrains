#!/bin/bash


#SBATCH --job-name=MFD_SMK          # Job name

#SBATCH --output=job_output_%j.log # Standard output and error log

#SBATCH --partition=cpu(all)       # Partition (queue) name

#SBATCH --nodes=1                  # Request two nodes

#SBATCH --ntasks-per-node=1        # Request 3 tasks per node

#SBATCH --cpus-per-task=4

#SBATCH --mem=16GB

# Load any necessary modules or environment settings here


# Your job commands go below this line

snakemake -s findviralstrains.smk --configfile config_files/no_ref_test.yml --cores 2 --rerun-incomplete
