import os
import re
import csv
import sys
import warnings
import resource
from datetime import datetime
from glob import glob

start_time = datetime.now()

print("  ______ _           ___      ___           _  _____ _             _            __      __           _                ___  __   ___  ")
print(" |  ____(_)         | \\ \\    / (_)         | |/ ____| |           (_)           \\ \\    / /          (_)              / _ \\/_ | / _ \\ ")
print(" | |__   _ _ __   __| |\\ \\  / / _ _ __ __ _| | (___ | |_ _ __ __ _ _ _ __  ___   \\ \\  / /__ _ __ ___ _  ___  _ __   | | | || || | | |")
print(" |  __| | | '_ \\ / _` | \\ \\/ / | | '__/ _` | |\\___ \\| __| '__/ _` | | '_ \\/ __|   \\ \\/ / _ \\ '__/ __| |/ _ \\| '_ \\  | | | || || | | |")
print(" | |    | | | | | (_| |  \\  /  | | | | (_| | |____) | |_| | | (_| | | | | \\__ \\    \\  /  __/ |  \\__ \\ | (_) | | | | | |_| || || |_| |")
print(" |_|    |_|_| |_|\\__,_|   \\/   |_|_|  \\__,_|_|_____/ \\__|_|  \\__,_|_|_| |_|___/     \\/ \\___|_|  |___/_|\\___/|_| |_|  \\___(_)_(_)___/ ")

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
###############
##   SETUP   ##
###############

# Get the current ulimit for file descriptors #
soft_limit, hard_limit = resource.getrlimit(resource.RLIMIT_NOFILE)

# Check if the soft limit is less than 2048 #
if soft_limit < 2048:
    warnings.warn("The current ulimit is less than 2048. This can cause Cuttlefish to fail.", ResourceWarning)
else:
    print("Ulimit is sufficient for program execution.")

# Check that the input folder exists
# Error if not
if os.path.isdir(READ_DIR) == False:
    raise OSError("\nRead directory " + READ_DIR + " not found\n")


samples = []
fastq_fullpath = []
fastq_filenames = []

for root, dirs, files in os.walk(READ_DIR):
    for name in files:
        # Remove the last 13 characters from each filename, will make it check for the .fastq ending later #
        sample_name = name[:-13]
        samples.append(sample_name)
        fastq_fullpath.append(os.path.join(root, name))
        fastq_filenames.append(sample_name)

if len(fastq_fullpath) < 1:
    raise OSError("\nNo .fastq files matching Illumina naming scheme found in input dir:\n" + READ_DIR + "\n")

# Get unique sample IDs
samples = list(set(samples))
samples.sort()

# Get expected number of lane files
# based on the sequencer type
def files_per_sequencer(sequencer):
    if (sequencer == "MiSeq"):
        return 1
    elif (sequencer == "NextSeq"):
        return 4
    else:
        raise NameError("\n" +sequencer + " is not a valid option for sequencer. Valid options:\nMiSeq, NextSeq\n")

expected_files = files_per_sequencer(SEQUENCER)

# A short function to add the output directory in front
# of filepaths, for ease of listing filenames below
def bd(filepath):
    return os.path.normpath(os.path.join(OUTPUT_DIR, ANALYSIS, filepath))

# Create the cuttlefish prefix
CF_PREF = bd("out")

fastq_filenames = set(fastq_filenames) # Deletes duplicate file entrys by converting to a set #
fastq_filenames = list(fastq_filenames)

fastq_filenames = [entry for entry in fastq_filenames if entry != ""] # Remake list with only populated values #
print(fastq_filenames)

######################
## HELPER FUNCTIONS ##
######################

# Functions for finding the fastq read files for a given sample
# First, a single function that finds the R1 and R2 files for a given sample
# and performs some checks to make sure they are properly paired
# and thebdre are the correct numbers
def find_read_files(wildcards):
    # Create search strings to match R1 and R2 files with standard Illumina names
    search_string_R1 = "^" + str({wildcards.sample}).strip("{'}") + "_.*_R1_.*\.fastq\.gz"
    search_string_R2 = "^" + str({wildcards.sample}).strip("{'}") + "_.*_R2_.*\.fastq\.gz"
    # Find the R1 and R2 files
    result_R1 = [f for f in fastq_fullpath if re.search(search_string_R1, os.path.basename(f))]
    result_R2 = [f for f in fastq_fullpath if re.search(search_string_R2, os.path.basename(f))]
    # Sort them
    result_R1.sort()
    result_R2.sort()

    # Check that R1 and R2 have an equal number of files
    # Give error if not (otherwise there will be issues with read pairing)
    if len(result_R1) != len(result_R2):
        raise OSError("\nUnequal number of R1 and R2 files found for sample: " + str({wildcards.sample}) + "\nR1 files: " + str(len(result_R1)) + "\nR2 files: " + str(len(result_R2)) + "\n")

    # Check that the the filenames for R1 and R2 match up properly
    # Throw error if not (out-of-order files will cause issues with read pairing)
    R1_check = [x.replace("_R1_", "") for x in result_R1]
    R2_check = [x.replace("_R2_", "") for x in result_R2]

    if R1_check != R2_check:
        print(result_R1)
        print(result_R2)
        raise OSError("\nFilenames of R1 and R2 reads for sample " + str({wildcards.sample}) + " do not match, this may cause issue with read pairing\nCheck above for the list of filenames causing this issue\n")

    # Check for the expected number of files
    # Print warning if not, but analysis can proceed
    if len(result_R1) != expected_files:
        print("Found " + str(len(result_R1)) + " sets of read files for sample " + str({wildcards.sample}).strip("{'}") + ", not " + str(expected_files))

    # Return lists of R1 and R2 files
    return [result_R1, result_R2]

# Then, two simple functions for use as input
# Which just return the R1 and R2 files, respectively
def find_R1_files(wildcards):
    R1_files = find_read_files(wildcards)[0]
    return R1_files
def find_R2_files(wildcards):
    R2_files = find_read_files(wildcards)[1]
    return R2_files

####################
## PIPELINE START ##
####################
#localrules: all

# Onstart and onsuccess, make a log files and copy the snakefile that was used
# into the results directory for posterity
#  TODO: add logging of all config file stuff??
# Add avoid over-writing, so that a pipeline can be run multiple times??
# Add in number of samples analyzed
pipeline_log_file = bd("pipeline_log.tsv")
onstart:
    os.makedirs(os.path.dirname(pipeline_log_file), exist_ok=True)
    with open(pipeline_log_file, 'w') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(["Raw data folder used:", READ_DIR])
        writer.writerow(["Prune count used:", PRUNE_COUNT])
        writer.writerow(["Start time:", start_time.strftime("%B %d, %Y: %H:%M:%S")])


onsuccess:
    smk_copy_command = 'cp findviralstrains.smk ' + str(bd("snakefile_used_copy.smk"))
    end_time = datetime.now()
    elapsed_time = end_time - start_time
    os.popen(smk_copy_command)
    with open(pipeline_log_file, 'a') as tsvfile:
        writer = csv.writer(tsvfile, delimiter='\t')
        writer.writerow(["End time:", end_time.strftime("%B %d, %Y: %H:%M:%S")])
        writer.writerow(["Elapsed time:", str(elapsed_time)])
        writer.writerow(["Copy of snakefile used stored at:", str(bd("snakefile_used_copy.smk"))])

# One rule to rule them all #
rule all:
	input:
		[bd(x) for x in expand("output_genomes/{input_list}/{input_list}_1_of_{numpaths}_vs_ref.txt", input_list=fastq_filenames, numpaths=["1", "2", "3"])]

rule trim_and_merge_raw_reads:
	input:
		raw_r1 = os.path.join(READ_DIR, "{sample}_R1_001.fastq"),
		raw_r2 = os.path.join(READ_DIR, "{sample}_R2_001.fastq"),
	output:
		trim_merged= (bd("processed_reads/trimmed/{sample}.merged.fq.gz")),
		trim_r1_pair= (bd("processed_reads/trimmed/{sample}.nomerge.pair.R1.fq.gz")),
		trim_r2_pair= (bd("processed_reads/trimmed/{sample}.nomerge.pair.R2.fq.gz")),
		trim_r1_nopair= (bd("processed_reads/trimmed/{sample}.nopair.R1.fq.gz")),
		trim_r2_nopair= (bd("processed_reads/trimmed/{sample}.nopair.R2.fq.gz")),
		rep_html= bd("logs/fastp/{sample}_trim_fastp.html"),
		rep_json= bd("logs/fastp/{sample}_trim_fastp.json")
#	threads: trim_threads # Unsure if this is needed #
	shell:
		"""
		fastp -i {input.raw_r1} -I {input.raw_r2} -m --merged_out {output.trim_merged} --out1 {output.trim_r1_pair} --out2 {output.trim_r2_pair} --unpaired1 {output.trim_r1_nopair} --unpaired2 {output.trim_r2_nopair} --detect_adapter_for_pe --cut_front --cut_front_window_size 5 --cut_front_mean_quality 20 -l 25 -j {output.rep_json} -h {output.rep_html} -w 1 2
		"""

# Unzip fastq files
rule Unzip:
    input:
        trim_merged = bd("processed_reads/trimmed/{sample}.merged.fq.gz"),
    output:
        unzipped = bd("processed_reads/trimmed/{sample}/{sample}.merged.fq"),
    shell:
        "gunzip -c {input.trim_merged} > {output.unzipped}"

# Make graph using BWT and our fm-index program
rule Create_graph:
    input:
        unzipped = bd("processed_reads/trimmed/{sample}/{sample}.merged.fq"),
    output:
        dbg = bd("dbg/{sample}/out.dbg"),
    params:
        pairdir = bd("processed_reads/trimmed/{sample}/"),
    shell:
        "target/release/assembly_graph_generator --input-dir {params.pairdir} --output-path {output.dbg} --kmer-len 27"

#Prune edges with small counts
rule Prune:
	input:
		dbg = bd("dbg/{sample}/out.dbg"),
	output:
		pruned_dbg = bd("dbg/{sample}/pruned/out.dbg"),
	shell:
		"python3 libs/prune/filter_reads.py {input.dbg} {output.pruned_dbg} {PRUNE_COUNT}"

rule Create_subgraphs:
    input:
        dbg = bd("dbg/{sample}/pruned/out.dbg"),
    output:
        graph_0 = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.dbg"),
        sources = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.sources"),
        sinks = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.sinks"),
        stats = bd("dbg/{sample}/out.dbg_subgraphs/graph_stats.txt"),
    shell:
        "target/release/graph_analyzer --dbg-file-name {input.dbg} --stats-output-file {output.stats}"

# Compress nodes with only one input and one output edge #
rule Compress:
	input:
		dbg = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.dbg"),
	output:
		comp_dbg = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0_compressed.dbg"),
	shell:
		"python3 libs/compress/compress.py {input.dbg} {output.comp_dbg}"

# Add super source and sink for ILP solver #
rule Add_super:
	input:
		comp_dbg = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0_compressed.dbg"),
		sources = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.sources"),
		sinks = bd("dbg/{sample}/pruned/out.dbg_subgraphs/graph_0.sinks"),
	output:
		swg = bd("wgs/super/{sample}.super.wg"),
	shell:
		"target/release/super_source_and_sink {input.sources} {input.sinks} {input.comp_dbg} {output.swg} graph_0"

# Uses Gurobi to try and sift our samples into different groups based on their reads #
rule Decompose:
	input:
		script = "libs/decompose/fracdecomp.py",
		swg = bd("wgs/super/{sample}.super.wg"),
	output:
		decomp = bd("decomp_results/{sample}.txt"),
		flow = bd("decomp_results/{sample}_1.paths"),
		flow2 = bd("decomp_results/{sample}_2.paths"),
		flow3 = bd("decomp_results/{sample}_3.paths"),
	shell:
		"python3 {input.script} -i {input.swg} -o {output.decomp} -M 3 --timelimit {DECOMP_TIME_LIMIT} -t {GUROBI_THREADS}"

# Runs rebuild.py to create a genome that follows the paths from Gurobi #
rule Rebuild_1:
	input:
		script = "libs/rebuild/rebuild.py",
		flow = bd("decomp_results/{sample}_1.paths"),
		swg = bd("wgs/super/{sample}.super.wg"),
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
		swg = bd("wgs/super/{sample}.super.wg"),
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
		swg = bd("wgs/super/{sample}.super.wg"),
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
