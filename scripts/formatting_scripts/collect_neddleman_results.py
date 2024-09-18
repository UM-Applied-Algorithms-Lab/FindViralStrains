import os
import csv

def process_files(output_genomes_path, decomp_results_path, output_file_path):
    # Dictionary to store data from output_genomes
    genome_data = {}

    # Traverse through the output_genomes directory
    for subdir in os.listdir(output_genomes_path):
        subdir_path = os.path.join(output_genomes_path, subdir)
        if os.path.isdir(subdir_path):
            for filename in os.listdir(subdir_path):
                if filename.endswith("_vs_ref.txt"):
                    sample_name = filename.replace("_vs_ref.txt", "")
                    
                    # Open the file and read the 25th and 26th lines
                    file_path = os.path.join(subdir_path, filename)
                    try:
                        with open(file_path, 'r') as file:
                            lines = file.readlines()
                            if len(lines) >= 26:
                                similarity_line = lines[24].strip()  # 25th line
                                gaps_line = lines[25].strip()       # 26th line
                                genome_data[sample_name] = (similarity_line, gaps_line)
                    except FileNotFoundError as e:
                        print(f"Error: {e}")

    # Write the extracted data to the output file
    with open(output_file_path, 'w') as output_file:
        output_file.write("Sample_Name, Similarity_Line, Gaps_Line\n")
        for sample_name, (similarity_line, gaps_line) in genome_data.items():
            output_file.write(f"{sample_name}, {similarity_line}, {gaps_line}\n")

    print(f"Data has been written to {output_file_path}")

if __name__ == "__main__":
    output_genomes_path = 'output_genomes'  # Replace if needed
    decomp_results_path = 'decomp_results'  # Replace with the actual path if needed
    output_file_path = 'extracted_data.csv'  # Output file path
    csv_file_path = output_file_path

    process_files(output_genomes_path, decomp_results_path, output_file_path)
