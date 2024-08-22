# Python env setup
import os
import re
import csv
import sys
from datetime import datetime

start_time = datetime.now()

print("  ______ _           ___      ___           _  _____ _             _            __      __           _                ___  __   ___  ")
print(" |  ____(_)         | \\ \\    / (_)         | |/ ____| |           (_)           \\ \\    / /          (_)              / _ \\/_ | / _ \\ ")
print(" | |__   _ _ __   __| |\\ \\  / / _ _ __ __ _| | (___ | |_ _ __ __ _ _ _ __  ___   \\ \\  / /__ _ __ ___ _  ___  _ __   | | | || || | | |")
print(" |  __| | | '_ \\ / _` | \\ \\/ / | | '__/ _` | |\\___ \\| __| '__/ _` | | '_ \\/ __|   \\ \\/ / _ \\ '__/ __| |/ _ \\| '_ \\  | | | || || | | |")
print(" | |    | | | | | (_| |  \\  /  | | | | (_| | |____) | |_| | | (_| | | | | \\__ \\    \\  /  __/ |  \\__ \\ | (_) | | | | | |_| || || |_| |")
print(" |_|    |_|_| |_|\\__,_|   \\/   |_|_|  \\__,_|_|_____/ \\__|_|  \\__,_|_|_| |_|___/     \\/ \\___|_|  |___/_|\\___/|_| |_|  \\___(_)_(_)___/ ")


#################
##   GLOBALS   ##
################# 
# Main config settings
ANALYSIS = config["analysis_ID"]
CONSENSUS_FILE = config["consensus_file"]
READ_DIR = config["read_dir"]
OUTPUT_DIR = config["output_dir"]
REF = config["ref_genome"]
SEQUENCER = config["sequencer"]
CONTIG_FILE = config["contig_file"]
CF_FILE = config["cf_file"]
###############
##   SETUP   ##
############### 

# Check that the input folder exists
# Error if not
if os.path.isdir(READ_DIR) == False:
    raise OSError("\nRead directory " + READ_DIR + " not found\n")

if os.path.isdir(READ_DIR) == False:
    raise OSError("\nConsensus file " + CONSENSUS_FILE + " not found\n")

# Pull all sample files from the raw data folder
# get filenames of the individual fastq files
fastq_fullpath = []
samples = []
fastq_filenames = [] # including the illumina stuff but not .fastq
for root, dirs, files in os.walk(READ_DIR):
    for name in files:
        if re.search("_S\d+_L\d+_R[12]_001.fastq$", name):
            samples.append(re.sub("_S\d+_L\d+_R\d+_001.fastq$", "", name))
            fastq_fullpath.append(os.path.join(root, name))
            fastq_filenames.append(re.sub(".fastq$", "", name))

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
        writer.writerow(["Consensus File used:", CONSENSUS_FILE])
        writer.writerow(["Reference genome used:", REF])
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

# all: The rule that looks for the final desired output files to initiate running all rules to generate those files.
rule all:
	input:
		expand(bd("output_genomes/{filename}.fasta"), filename=fastq_filenames)

# Cleans files and makes sure to remove trailing whitespace that may cause issues later #
rule create_contigs: 
	input:
		file = CONSENSUS_FILE,
		script = "scripts/createcontigs.sh",
	output:
		file = bd("all_samples_consensus_contigs.fasta")
	shell:
		"bash {input.script} {input.file} {output.file}"	

# Runs Cuttlefish #
rule Cuttlefish: 
	input:
		file = bd("all_samples_consensus_contigs.fasta")
	output:
		seg = bd("out.cf_seg"), 
		seq = bd("out.cf_seq") 
	shell:
		#"ulimit -n 2048 ",
		"rm -f " + bd("out.json") + " && cuttlefish build -s {input.file} -t 1 -o {CF_PREF} -f 3 -m 12"

rule Mer_graph: 
	input:
		script = "scripts/edgemer.py",
		seg = bd("out.cf_seg"), 
		seq = bd("out.cf_seq"), 
	output:
		file = bd("out.mg"),
	shell:
		"python3 {input.script} -k 27 -c {CF_PREF} -o {output.file}"

rule Run_jf:
    input:
        script = "scripts/runjf.sh",
        mg = bd("out.mg"),
        reads = READ_DIR + "{sample}.fastq"
    output:
        bd("wgs/{sample}.wg")
    shell:
        "{input.script} {input.reads} {input.mg} {output}"


# Uses Gurobi to try and sift our samples into different groups based on their reads #
rule Decompose:
	input:
		wg = bd("wgs/{sample}.wg"),
		script = "scripts/fracdecomp.py",
	output:
		decomp = bd("decomp_results/{sample}.txt")
	shell:
		"python3 {input.script} -i {input.wg} -o {output.decomp} -M 3" # TODO Change name scheme #

# Future rule to be added to use format_to_graph that will create graphs showing each path #

# Runs rebuild.sh to create a genome that follows each of the paths #
rule Rebuild:
	input:
		script = ("scripts/rebuild.sh"),
		flow = bd("decomp_results/{sample}.txt_1_paths.txt"),
		cf_seg = bd("out.cf_seg"),
	output:
		genome = (bd("output_genomes/{sample}.fasta"))
	shell:
		"bash {input.script} {input.flow} {input.cf_seg} {output.genome}"


