#!/bin/bash

# Input and output directories
INPUT_DIR="/home/mikhail/TableDir"
FINAL_OUTPUT_FILE="${INPUT_DIR}/latex_table.tex"

# Clear the final output file if it exists
> "$FINAL_OUTPUT_FILE"

# Start the LaTeX table
echo "\\begin{tabular}{|c|c|c|c|c|c|c|c|c|c|}" >> "$FINAL_OUTPUT_FILE"
echo "\\hline" >> "$FINAL_OUTPUT_FILE"

# Add table headers
echo "Sample & Error-1path & Path1-wt & Error-2paths & Path1-wt & Path2-wt & Error-3paths & Path1-wt & Path2-wt & Path3-wt \\\\" >> "$FINAL_OUTPUT_FILE"
echo "\\hline" >> "$FINAL_OUTPUT_FILE"

# Counter for sample numbers
counter=1

# Loop over each file ending in R1_001.txt in the input directory
for file in "$INPUT_DIR"/*R1_001.txt; do
    if [ -f "$file" ]; then
        # Read the single line from the file
        line=$(cat "$file")
        
        # Replace the first item with "Sample N" and format numbers
        formatted_line="Sample $counter & $(echo "$line" | awk '{ for (i=2; i<=11; i++) printf "%.2f ", $i }' | sed 's/ / \& /g')"
        
        # Ensure the line has exactly ten columns by truncating if necessary
        formatted_line=$(echo "$formatted_line" | cut -d'&' -f1-10)
        
        # Append the formatted line to the output file
        echo "$formatted_line \\\\" >> "$FINAL_OUTPUT_FILE"
        echo "\\hline" >> "$FINAL_OUTPUT_FILE"
        
        # Increment the counter
        ((counter++))
    fi
done

# End the LaTeX table
echo "\\end{tabular}" >> "$FINAL_OUTPUT_FILE"

echo "LaTeX table has been written to $FINAL_OUTPUT_FILE"

