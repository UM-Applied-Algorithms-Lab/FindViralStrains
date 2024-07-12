#!/bin/bash
file="$1"
InputMergraph="$2"
OutputFile="$3"
start_time=$(date +%s) 

# Function to calculate and print the elapsed time
print_elapsed_time() {
    local current_time=$(date +%s)
    local elapsed_time=$((current_time - start_time))
    echo "Elapsed time: ${elapsed_time} seconds"
}

# files are called SAMPLE_NAME.fastq
echo "Processing file: $file"
print_elapsed_time
JellyfishIndex="$OutputFile.jf"
jellyfish count -m 28 -s 100M -t 1 -C $file -o $JellyfishIndex

echo "jellyfish index built for $file"
print_elapsed_time

# Open and read the input file line by line 
while IFS= read -r line; do
	# Remove leading/trailing whitespace
	line=$(echo "$line" | xargs)

	# Split the line into an array based on whitespace
	read -a array <<< "$line"
		
	# Get the number of fields in the array
	fields=${#array[@]}
	
	# if it's a line for an edge
	if [ "$fields" -eq 3 ]; then
		seq="${array[2]}"
		# query the index for the edgemer count
		count=$(jellyfish query $JellyfishIndex "$seq" | awk '{print $2}')
		# print
		echo -e "${array[0]}\t${array[1]}\t$count" >> $OutputFile
	else
		echo $line >> $OutputFile
	fi

done < "$InputMergraph"

echo "wg created for $file"
print_elapsed_time
