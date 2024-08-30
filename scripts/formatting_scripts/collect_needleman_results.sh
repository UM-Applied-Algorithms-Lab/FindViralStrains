#!/bin/bash

# Check if both directories are provided as arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <needle_dir> <paths_dir>"
    exit 1
fi

# Starting directories
needle_dir="$1"
paths_dir="$2"

# Output file to store results
output_file="similarity_gap_results.txt"
> "$output_file" # Clear the file if it exists

# Temporary file for storing results before appending path lines
temp_output_file="temp_similarity_gap_results.txt"
> "$temp_output_file"

# Traverse through subdirectories and process each _vs_ref.txt file
find "$needle_dir" -type f -name "*_vs_ref.txt" | while read -r file; do
    echo "Processing file: $file"

    # Extract Similarity and Gaps information from the file
    similarity=$(grep -oP 'Similarity:\s*\d+/\d+\s*\(\K[0-9.]+(?=%\))' "$file")
    gaps=$(grep -oP 'Gaps:\s*\d+/\d+\s*\(\K[0-9.]+(?=%\))' "$file")

    # Get the base name without _vs_ref.txt
    base_name=$(basename "$file" "_vs_ref.txt")

    # Write the results to the temporary output file
    echo "$base_name: Similarity=${similarity}%, Gaps=${gaps}%" >> "$temp_output_file"
done

# Parse the second directory for .paths files and append the Objective Value
find "$paths_dir" -type f -name "*1.paths" | while read -r paths_file; do
    echo "Processing paths file: $paths_file"

    # Extract the Objective Value using grep
    objective_value=$(grep -oP '^Objective Value:\s*\K[0-9.]+$' "$paths_file")

    # Get the base name to match with the first part and remove trailing underscore
    paths_base_name=$(basename "$paths_file" "1.paths")
    paths_base_name="${paths_base_name%_}" # Remove trailing underscore

    # Append the Objective Value to the corresponding entry in the temp output file
    sed -i "s/^$paths_base_name: /&Objective_Value=${objective_value}, /" "$temp_output_file"
done

# Sort the output file by the Objective Value (lowest to highest)
#awk -F'Objective_Value=' '{print $2, $0}' "$temp_output_file" | sort -n | cut -d' ' -f2- | uniq > "$output_file"

echo "Results have been written to $output_file"

