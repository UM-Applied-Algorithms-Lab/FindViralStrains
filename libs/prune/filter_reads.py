import sys

# Input and output file paths from command-line arguments
dbg_file = sys.argv[1]
filtered_dbg_file = sys.argv[2]  # Output file for filtered mg data
#num_prune = sys.argv[3] --temp for testing purposes

# Process the mg_file to filter based on the count (third number)
with open(dbg_file, 'r') as file:
    with open(filtered_dbg_file, 'w') as outfile:
        for line in file:
            parts = line.split()  # Split the line into parts
            
            # Check if the line has at least 3 parts
            if len(parts) >= 4:
                first_argument = parts[0]
                second_argument = parts[1]
                count_argument = parts[2]
                
                try:
                    # Try to convert the count argument to a float
                    count_value = float(count_argument)
                    
                    # Skip this line if the count value is less than 25
                    if count_value < 15:
                        continue

                except ValueError:
                    # If it's not a number, continue with the line as is
                    pass

            # Write the line to the output file if it passed the check
            outfile.write(line)
