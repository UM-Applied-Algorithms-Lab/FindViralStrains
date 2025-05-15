#!/bin/bash
apptainer run ../findviralstrains.sif

for sample_dir in test/data/NoRefTest/dbg/*/; do
    sample=$(basename "$sample_dir")
    input_file="test/data/NoRefTest/dbg/$sample/out.dbg"
    job_script="batch_scripts/prune_job_${sample}.sh"

    sed \
        -e "s|{{sample}}|$sample|g" \
        -e "s|{{input_file}}|$input_file|g" \
	batchloop.sh | sbatch
done
