#!/bin/bash
Inputfile="$1"
Outputfile="$2"

rm -f $Outputfile




# Open and read the input file line by line 
while IFS= read -r line; do
	# Check if line starts with a '>' symbol (then it's a name line_
	if [[ $line == \>* ]]; then
		name=$line
	else
		# otherwise, it's a sequence line
		# Remove leading/trailing whitespace and replace Ns with spaces
		line=$(echo "$line" | xargs)
	    echo $line | sed 's/N/ /g' | awk -v header="$name" '{for (i=1; i<=NF; i++) if ($i != "") {print header "_" i; print $i}}' >> $Outputfile
	 fi
done < "$Inputfile"


