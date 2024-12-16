import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import sys
import glob

output_file = sys.argv[1] if len(sys.argv) > 1 else None

# List to store the third arguments
graph_weights = []

# Loop through all `.wg` files in the current directory
for input_file in glob.glob("*.wg"):
    with open(input_file, 'r') as file:
        for line in file:
            parts = line.split()
            if len(parts) >= 3:
                try:
                    graph_weights.append(float(parts[2]))
                except ValueError:
                    continue

# Check if data is available
if not graph_weights:
    print("No valid data found in .wg files.")
    sys.exit(1)

# Determine the range of the data
min_value = min(graph_weights)
max_value = max(graph_weights)

# Create bins with a step of 25
bin_edges = list(range(int(min_value), int(max_value) + 25, 25))

# Plot the histogram
plt.figure(figsize=(12, 6))  # Increased figure size for readability
plt.hist(graph_weights, bins=bin_edges, edgecolor='black', color='skyblue', align='mid')

# Add titles and labels with detailed information
plt.title('Histogram of Weights on Graph Edges', fontsize=16)
plt.xlabel('Weight Amount', fontsize=14)
plt.ylabel('Number of Occurrences', fontsize=14)

# Apply a logarithmic scale to the Y-axis
plt.yscale('log')

# Customize the X-axis ticks to match the bins
plt.gca().xaxis.set_major_locator(ticker.MultipleLocator(100))  # Set major ticks every 100 units
plt.grid(axis='y', linestyle='--', alpha=0.7)

# Rotate X-axis labels vertically
plt.xticks(rotation=90, fontsize=7)
plt.yticks(fontsize=12)

# Display or save the plot
plt.tight_layout()

if output_file:
    plt.savefig(output_file)
    print(f"Graph saved to {output_file}")
else:
    plt.show()
