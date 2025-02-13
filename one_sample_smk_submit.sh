#!/bin/bash


#SBATCH --job-name=MFD_SMK          # Job name

#SBATCH --output=slurm_out/job_output_%j.log # Standard output and error log
#SBATCH --error=slurm_out/job_output_%j.err  # Standard output and error log

#SBATCH --partition=cpu(all)       # Partition (queue) name

#SBATCH --nodes=1                  # Request two nodes

#SBATCH --ntasks-per-node=1        # Request 3 tasks per node

#SBATCH --cpus-per-task=4

#SBATCH --mem=16GB

# Load any necessary modules or environment settings here


# Your job commands go below this line

ulimit -n 2048

## Enable Conda environment
source /hellgate/home/$USER/.bashrc
conda activate Test

## Run pipeline
snakemake -s findviralstrains.smk --configfile config_files/one_sample.yml --cores 2 --rerun-incomplete
