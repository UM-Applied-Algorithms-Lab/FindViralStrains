#!/bin/bash

# File paths (adjust as necessary)

#"~/Code/MFD-ILP/FindViralStrains/output/NoRefTest/out.cf_seg"
path_file="$1"
seq_file="$2"
# Add seg file and ref genome with bd() and output #
ref_genome="$3"
#"~/Code/MFD-ILP/FindViralStrains/reference_genomes/covid19ref.fasta"
bd_outfile="$4"
# Extract the base name without the extension and directory
base_name=$(basename "$path_file" | cut -d'.' -f1)

# Construct the output directory with the part of the name and .fasta
output_dir="${bd_outfile}/output_genomes/${base_name}"

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Function to get the reverse complement of an RNA sequence
reverse_complement() {
  local seq="$1"
  # Reverse the sequence and complement it
  echo "$seq" | rev | tr 'ACGTacgt' 'TGCAtgca'
}

# Step 1: Read sequences into an associative array
declare -A sequences
while IFS=$'\t' read -r node_id sequence; do
  sequences["$node_id"]="$sequence"
done < "$seq_file"

# Step 2: Count total number of alignments
total_alignments=$(grep -c '^[0-9]\+' "$path_file")

# Step 3: Process the paths and reconstruct genomes
counter=1
while IFS= read -r line; do
  if [[ "$line" =~ ^[0-9]+(\.[0-9]+)? ]]; then
    # This is a line with a weight and a path
    weight=$(echo "$line" | cut -d' ' -f1)
    path=$(echo "$line" | cut -d' ' -f2-)

    # Initialize final genome
    genome=""
    is_first_node=true

    # Process each node in the path
    IFS=' ' read -ra nodes <<< "$path"
    for node in "${nodes[@]}"; do
      # Check the sign and clean the node identifier
      if [[ "$node" == *"-" ]]; then
        clean_node=$(echo "$node" | sed 's/-//g')
        sequence="${sequences[$clean_node]}"
        if [[ -n "$sequence" ]]; then
          # Handle the overlap by deleting the first 27 characters for non-first nodes
          if [ "$is_first_node" = false ]; then
            sequence="${sequence:27}"
          fi
          genome+=$(reverse_complement "$sequence")
        else
          echo "Warning: Node $clean_node not found in sequences."
        fi
      else
        clean_node=$(echo "$node" | sed 's/+//g')
        sequence="${sequences[$clean_node]}"
        if [[ -n "$sequence" ]]; then
          # Handle the overlap by deleting the first 27 characters for non-first nodes
          if [ "$is_first_node" = false ]; then
            sequence="${sequence:27}"
          fi
          genome+="$sequence"
        else
          echo "Warning: Node $clean_node not found in sequences."
        fi
      fi
      is_first_node=false
    done

    # Write to a file in the output directory
    output_file="$output_dir/genome_${counter}_of_${total_alignments}.fasta"
    echo -e ">Weight: $weight\n$genome" > "$output_file"

    # Perform alignment with the reference genome
    alignment_file="$output_dir/alignment_${counter}_of_${total_alignments}.txt"
#    needle -asequence "$ref_genome" -bsequence "$output_file" -gapopen 10 -gapextend 0.5 -outfile "$alignment_file"
# TODO Add the step that will do this
    ((counter++))
  fi
done < "$path_file"

