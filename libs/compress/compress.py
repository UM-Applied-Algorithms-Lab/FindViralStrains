with open('/home/mikhailu/Projects/FindViralStrains/output/FakeDataTest/dbg/testdata_L001/out.dbg', 'r') as file:
    lines = file.readlines()

# Split each line into columns
data = []
for line in lines:
    columns = line.strip().split('\t')
    if columns:  # skip empty lines
        data.append((int(columns[0]), '\t'.join(columns[1:])))  # store as tuple (first_col, rest_of_line)

# Sort by the first column
sorted_data = sorted(data, key=lambda x: x[0])

# Write sorted 
with open('sorted_output.txt', 'w') as file:
    for item in sorted_data:
        file.write(f"{item[0]}\t{item[1]}\n")

# Adjacency Lists and checking for size of those lists
# Make if exception for weight 0 edges (super source and sink) to not be compressed