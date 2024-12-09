#!/bin/bash

# Path to the input directory
INPUT_DIR="../../../BaseCalls"
# Path to the output directory
OUTPUT_DIR="./aligned_outputs"
# HISAT2 index prefix
INDEX="covid_index"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Loop through R1 files in the directory
for R1_FILE in "$INPUT_DIR"/*R1_001.fastq; do
    # Derive the corresponding R2 file
    R2_FILE="${R1_FILE/_R1_/_R2_}"

    # Extract the sample name (e.g., E2081_S71_L001 from E2081_S71_L001_R1_001.fastq)
    SAMPLE_NAME=$(basename "$R1_FILE" | sed 's/_R1_001.fastq//')

    # Define output file name
    OUTPUT_FILE="$OUTPUT_DIR/${SAMPLE_NAME}_aligned.sam"

    # Run HISAT2 for the sample with your specified parameters
    hisat2 -x "$INDEX" \
        -1 "$R1_FILE" \
        -2 "$R2_FILE" \
        -S "$OUTPUT_FILE" \
        --threads 8 \
        --score-min L,0,-0.2 \
        --pen-canintronlen G,-8,1

    echo "Finished processing $SAMPLE_NAME"
done

