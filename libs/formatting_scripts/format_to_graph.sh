#!/bin/bash

# Input file
input_file="input.txt"

# Output DOT file
dot_file="graph.dot"

# Output PDF file
pdf_file="graph.pdf"

# Initialize the DOT file
echo "digraph G {" > $dot_file

# Read the input file and process each line
while read -r line; do
    # Extract the first two columns (nodes)
    node1=$(echo $line | awk '{print $1}')
    node2=$(echo $line | awk '{print $2}')
    
    # Add the edge to the DOT file
    echo "    $node1 -> $node2;" >> $dot_file
done < $input_file

# Close the DOT file
echo "}" >> $dot_file

# Generate the PDF using dot
dot -Tpdf $dot_file -o $pdf_file

echo "Graph generated as $pdf_file"