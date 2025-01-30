import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys

# Input and output file paths from command-line arguments
wg_file = sys.argv[1]
mg_file = sys.argv[2]
output_file = sys.argv[3]  # Output file for filtered wg data
filtered_mg_file = sys.argv[4]  # Output file for filtered mg data

# Store valid points from wg_file
valid_points = set()

# Process the wg_file
with open(wg_file, 'r') as file:
    with open(output_file, 'w') as outfile:
        for line in file:
            parts = line.split()  # Split the line into parts
            
            # Check if the line has at least 3 parts
            if len(parts) >= 3:
                first_argument = parts[0]
                second_argument = parts[1]
                third_argument = parts[2]
                
                try:
                    # Try to convert the third argument to a float
                    third_value = float(third_argument)
                    
                    # Skip this line if the third value is less than 20
                    if third_value < 25:
                        continue

                    # Add valid (first, second) pair to the set
                    valid_points.add((first_argument, second_argument))
                except ValueError:
                    # If it's not a number, continue with the line as is
                    pass

            # Write the line to the output file if it passed the check
            outfile.write(line)

# Process the mg_file to match the valid points from wg_file
with open(mg_file, 'r') as file:
    with open(filtered_mg_file, 'w') as outfile:
        for line in file:
            parts = line.split()  # Split the line into parts
            
            # Check if the line has at least 2 parts
            if len(parts) >= 2:
                first_argument = parts[0]
                second_argument = parts[1]

                # Check if the (first, second) pair is in the valid points set
                if (first_argument, second_argument) in valid_points:
                    outfile.write(line)
