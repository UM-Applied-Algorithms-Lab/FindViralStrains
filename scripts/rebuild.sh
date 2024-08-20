#!/bin/bash

# File paths (adjust as necessary) #
seq_file="out.cf_seg"
path_file="decomp_results/54754002_S19_L001_R2_001.txt_1_paths.txt"
output_dir="output_genomes"

# Create output directory if it doesn't exist #
mkdir -p "$output_dir"

# Function to get the reverse complement of a RNA sequence #
reverse_complement() {
    local seq="$1"
    # Reverse the sequence and complement it #
    echo "$seq" | rev | tr 'ACGTacgt' 'TGCAtgca'
}

# Step 1: Read sequences into an associative array #
declare -A sequences
while IFS=$'\t' read -r node_id sequence; do
    sequences["$node_id"]="$sequence"
done < "$seq_file"

# Step 2: Process the paths and reconstruct genomes #
counter=1
while IFS= read -r line; do
    if [[ "$line" =~ ^[0-9]+(\.[0-9]+)? ]]; then
        # This is a line with a weight and a path
        weight=$(echo "$line" | cut -d' ' -f1)
        path=$(echo "$line" | cut -d' ' -f2-)

        # Initialize final genome #
        genome=""
        first_node=true

        # Process each node in the path #
        IFS=' ' read -ra nodes <<< "$path"
        for node in "${nodes[@]}"; do
            # Check the sign and clean the node identifier #
            if [[ "$node" == *"-" ]]; then
                clean_node=$(echo "$node" | sed 's/-//g')
                if [[ -n "${sequences[$clean_node]}" ]]; then
                    # Add the reverse complement if needed #
                    sequence=$(reverse_complement "${sequences[$clean_node]}")
                else
                    echo "Warning: Node $clean_node not found in sequences."
                    continue
                fi
            else
                clean_node=$(echo "$node" | sed 's/+//g')
                if [[ -n "${sequences[$clean_node]}" ]]; then
                    # Add the sequence directly #
                    sequence="${sequences[$clean_node]}"
                else
                    echo "Warning: Node $clean_node not found in sequences."
                    continue
                fi
            fi

            # Check if it's the first node
            if $first_node; then
                genome+="$sequence"
                first_node=false
            else
                # Cut off the first 27 characters for subsequent nodes
                genome+="${sequence:27}"
            fi
        done

        # Write to a file in the output directory #
        output_file="$output_dir/genome_${counter}.fasta"
        echo -e ">Weight: $weight\n$genome" > "$output_file"
        ((counter++))
    fi
done < "$path_file"

