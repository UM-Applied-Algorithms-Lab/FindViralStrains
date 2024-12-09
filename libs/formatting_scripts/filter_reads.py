import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]  # Output file specified by the second argument

graph_weights = []

with open(input_file, 'r') as file:
    # Open the output file in write mode
    with open(output_file, 'w') as outfile:
        # Loop through each line in the input file
        for line in file:
            # Split the line into parts (assuming whitespace as the delimiter)
            parts = line.split()
            
            # Check if the line has at least 3 parts
            if len(parts) >= 3:
                third_argument = parts[2]
                
                try:
                    # Try to convert the third argument to a float
                    third_value = float(third_argument)
                    
                    # If the third value is less than 25, skip this line entirely
                    if third_value < 10:
                        continue  # Skip writing this line to the output file
                except ValueError:
                    # If it's not a number, continue with the line as is
                    pass
            
            # If the line passed the check, write it to the output file
            outfile.write(line)
