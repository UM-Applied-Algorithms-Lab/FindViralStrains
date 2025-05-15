#!/bin/bash
#SBATCH --job-name=prune_{{sample}}
#SBATCH --output=logs/prune_{{sample}}.out
#SBATCH --error=logs/prune_{{sample}}.err
#SBATCH --time=03:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=2

sample="B0841_S78_L001"
input_file="test/data/NoRefTest/dbg/$sample/out.dbg"
for i in {0..30}; do
    pruned_dir="test/data/NoRefTest/dbg/$sample/pruned_$i"
    pruned_file="$pruned_dir/out.dbg"
    pruned_ref_file="$pruned_dir/ref_out.dbg" 
    mkdir -p "$pruned_dir"

    python3 libs/prune/filter_reads.py "$input_file" "$pruned_file" "$i"
    python3 libs/prune/add_reference.py "$pruned_file" "/mnt/beegfs/projects/tb208541/data_from_lucy/reference_genome/GCF_009858895.2_ASM985889v3_genomic_notail.fna" "$pruned_ref_file" 
    stats_dir="dbg_parallel/$sample"
    stats_file="$stats_dir/graph_stats_$i.txt"
    mkdir -p "$stats_dir"

    apptainer run --bind /mnt/beegfs/projects/tb208541  ../findviralstrains.sif  target/release/graph_analyzer --dbg-file-name "$pruned_ref_file" --node-percent-cutoff 0 --stats-output-file "$stats_file"
done

