import os
import re
import csv
import sys
import warnings
import resource
from datetime import datetime
from glob import glob
from pathlib import Path


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

def get_subgraph_indices(sample):
    """
    Reads the output directory from config and returns a list of all subgraph indices.
    """
    # Construct the path to the subgraphs directory using config OUTPUT_DIR
    subgraphs_dir = Path(OUTPUT_DIR) / "graphs" / sample / "pruned.dbg_subgraphs"
    
    # Find all .dbg files in the directory
    dbg_files = glob.glob(str(subgraphs_dir / "graph_*.dbg"))
    
    # Extract the numeric indices from the filenames
    subgraph_indices = []
    for file_path in dbg_files:
        try:
            # Extract the number between graph_ and .dbg
            base = os.path.basename(file_path)
            num = int(base.split('_')[1].split('.')[0])
            subgraph_indices.append(num)
        except (IndexError, ValueError):
            continue
    
    # Return sorted unique indices
    return sorted(subgraph_indices)


	# One rule to rule them all #
rule all:
	input:
		[bd(x) for x in expand("output_genomes/{input_list}/{input_list}_1_of_{numpaths}_vs_ref.txt", input_list=fastq_filenames, numpaths=["1", "2", "3"])]

	
# Compress nodes with only one input and one output edge #
rule Compress:
	input:
		dbg = bd("graphs/{sample}/pruned.dbg_subgraphs/graph_0.dbg"),
	output:
		comp_dbg = bd("graphs/{sample}/out.dbg_subgraphs/graph_0_compressed.dbg"),
	shell:
		"python3 libs/compress/compress.py {input.dbg} {output.comp_dbg}"

# Add super source and sink for ILP solver #
rule Add_super:
	input:
		comp_dbg = bd("graphs/{sample}/out.dbg_subgraphs/graph_0_compressed.dbg"),
		sources = bd("graphs/{sample}/pruned.dbg_subgraphs/graph_0.sources"),
		sinks = bd("graphs/{sample}/pruned.dbg_subgraphs/graph_0.sinks"),
	output:
		swg = bd("graphs/{sample}.super.dbg"),
	shell:
		"target/release/super_source_and_sink {input.sources} {input.sinks} {input.comp_dbg} {output.swg} graph_0"

# Uses Gurobi to try and sift our samples into different groups based on their reads #
rule Decompose:
	input:
		script = "libs/decompose/kleast_errors.py",
		swg = bd("graphs/{sample}.super.dbg"),
	output:
		flow = bd("decomp_results/{sample}_1.paths"),
		flow2 = bd("decomp_results/{sample}_2.paths"),
		flow3 = bd("decomp_results/{sample}_3.paths"),
	params:
		decomp = bd("decomp_results/{sample}.txt"),
	shell:
		"python3 {input.script} -i {input.swg} -o {params.decomp} -M 3 --timelimit {DECOMP_TIME_LIMIT} -t {GUROBI_THREADS} --visualize {VISUALIZE}"

# Runs rebuild.py to create a genome that follows the paths from Gurobi #
rule Rebuild_1:
	input:
		script = "libs/rebuild/rebuild.py",
		flow = bd("decomp_results/{sample}_1.paths"),
		swg = bd("graphs/{sample}.super.dbg"),
	output:
		genome = bd("output_genomes/{sample}/{sample}_1_of_1.fasta"),
	params:
		outtemp = bd("output_genomes/{sample}/{sample}.fasta")
	shell:
		"""
		python3 {input.script} {input.flow} {input.swg} {params.outtemp}
		"""

rule Rebuild_2:
	input:
		script = "libs/rebuild/rebuild.py",
		flow = bd("decomp_results/{sample}_2.paths"),
		swg = bd("graphs/{sample}.super.dbg"),
	output:
		genome = bd("output_genomes/{sample}/{sample}_1_of_2.fasta"),
		genome2 = bd("output_genomes/{sample}/{sample}_2_of_2.fasta"),
	params:
		outtemp = bd("output_genomes/{sample}/{sample}.fasta")
	shell:
		"""
		python3 {input.script} {input.flow} {input.swg} {params.outtemp}
		"""

rule Rebuild_3:
	input:
		script = "libs/rebuild/rebuild.py",
		flow = bd("decomp_results/{sample}_3.paths"),
		swg = bd("graphs/{sample}.super.dbg"),
	output:
		genome = bd("output_genomes/{sample}/{sample}_1_of_3.fasta"),
		genome2 = bd("output_genomes/{sample}/{sample}_2_of_3.fasta"),
		genome3 = bd("output_genomes/{sample}/{sample}_3_of_3.fasta"),
	params:
		outtemp = bd("output_genomes/{sample}/{sample}.fasta")
	shell:
		"""
		python3 {input.script} {input.flow} {input.swg} {params.outtemp}
		"""

# Compares our newly constructed genomes to original covid reference using Needleman-Wunsch #
# A more modern reference could be used, or regional samples as well #
# There are also better versions of this, and I may make my own #
rule Compare_1:
	input:
		rebuilt_genome = bd("output_genomes/{sample}/{sample}_1_of_1.fasta"),
		origin_covid = ("reference_genomes/covid19ref.fasta")
	output:
		compar_file = bd("output_genomes/{sample}/{sample}_1_of_1_vs_ref.txt")
	shell:
		"needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file}"

# Compares genomes from the two path result to the reference #
rule Compare_2:
	input:
		rebuilt_genome_1 = bd("output_genomes/{sample}/{sample}_1_of_2.fasta"),
		rebuilt_genome_2 = bd("output_genomes/{sample}/{sample}_2_of_2.fasta"),
		origin_covid = ("reference_genomes/covid19ref.fasta")
	output:
		compar_file_1 = bd("output_genomes/{sample}/{sample}_1_of_2_vs_ref.txt"),
		compar_file_2 = bd("output_genomes/{sample}/{sample}_2_of_2_vs_ref.txt")
	shell:
		"""
		needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome_1} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file_1}
		needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome_2} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file_2}
		"""

# Compares genomes from the three path result to the reference #
rule Compare_3:
	input:
		rebuilt_genome_1 = bd("output_genomes/{sample}/{sample}_1_of_3.fasta"),
		rebuilt_genome_2 = bd("output_genomes/{sample}/{sample}_2_of_3.fasta"),
		rebuilt_genome_3 = bd("output_genomes/{sample}/{sample}_3_of_3.fasta"),
		origin_covid = ("reference_genomes/covid19ref.fasta")
	output:
		compar_file_1 = bd("output_genomes/{sample}/{sample}_1_of_3_vs_ref.txt"),
		compar_file_2 = bd("output_genomes/{sample}/{sample}_2_of_3_vs_ref.txt"),
		compar_file_3 = bd("output_genomes/{sample}/{sample}_3_of_3_vs_ref.txt")
	shell:
		"""
		needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome_1} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file_1}
		needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome_2} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file_2}
		needle -asequence {input.origin_covid} -bsequence {input.rebuilt_genome_3} -gapopen 10 -gapextend 0.5 -outfile {output.compar_file_3}
		"""
