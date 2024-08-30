#!/bin/bash

# Check if the directory is provided as an argument #
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <directory-path>"
    exit 1
fi

# Starting directory containing subdirectories #
root_dir="$1"

# Output file to store results #
output_file="similarity_gap_results.txt"
> "$output_file"  # Clear the file if it exists #

# Traverse through subdirectories and process each _vs_ref.txt file #
find "$root_dir" -type f -name "*_vs_ref.txt" | while read -r file; do
    echo "Processing file: $file"

    # Extract Similarity and Gaps information from the file #
    similarity=$(grep -oP 'Similarity:\s*\d+/\d+\s*\(\K[0-9.]+(?=%\))' "$file")
    gaps=$(grep -oP 'Gaps:\s*\d+/\d+\s*\(\K[0-9.]+(?=%\))' "$file")

    # Write the results to the output file #
    echo "$(basename "$file"): Similarity=${similarity}%, Gaps=${gaps}%" >> "$output_file"
done

echo "Results have been written to $output_file"

