import os
import re
import matplotlib.pyplot as plt

def process_file(file_path):
    gap_locations = ["/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/output_genomes"]
    in_alignment_section = False
    
    with open(file_path, 'r') as file:
        lines = file.readlines()
        for line in lines:
            if line.startswith("NC_045512.2"):
                in_alignment_section = True
                continue
            
            if line.startswith("Weight"):
                in_alignment_section = False
                continue
            
            if in_alignment_section:
                # Extract gap locations from alignment lines
                gap_matches = [m.start() for m in re.finditer(r'-', line)]
                gap_locations.extend(gap_matches)
    
    print(f"Gaps found in {file_path}: {gap_locations}")
    return gap_locations

def plot_gaps(all_gap_locations):
    plt.figure(figsize=(10, 6))
    plt.scatter(range(len(all_gap_locations)), all_gap_locations, c='blue', marker='o')
    plt.title('Gap Locations')
    plt.xlabel('Index')
    plt.ylabel('Position')
    plt.grid(True)
    plt.savefig('gaps_plot.png')
    plt.show()

def main(directory_path):
    all_gap_locations = []

    for filename in os.listdir(directory_path):
        if filename.endswith('_vs_ref.txt'):
            file_path = os.path.join(directory_path, filename)
            print(f"Processing file: {file_path}")
            gap_locations = process_file(file_path)
            all_gap_locations.extend(gap_locations)
    
    plot_gaps(all_gap_locations)

if __name__ == "__main__":
    directory_path = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/output_genomes"  # Replace with your directory path
    main(directory_path)
