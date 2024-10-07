import os
import csv
import fnmatch

excel_doc = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/extracted_data.csv"
paths_dir = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/decomp_results/"
updated_excel_doc = "/home/mikhail/Code/MFD-ILP/FindViralStrains/output/NoRefTest/extracted_data_updated.csv"

# Check if the Excel document exists
if os.path.exists(excel_doc):
    with open(excel_doc, mode='r', newline='', encoding='utf-8') as file:
        csv_reader = csv.reader(file)
        rows = list(csv_reader)  # Store the rows to modify

    with open(updated_excel_doc, mode='w', newline='', encoding='utf-8') as updated_file:
        csv_writer = csv.writer(updated_file)

        # Iterate through rows in the original CSV
        for row in rows:
            if row:  # Ensure the row is not empty
                first_word = row[0].split()[0]  # Split the first column
                sample_id = first_word[:-7]  # Remove the last seven characters
                path = first_word[-1]  # Save the last character

                # Print Sample ID and Path
                print(f"Sample ID: {sample_id}, Path: {path}")

                # Look for files in paths_dir matching the pattern
                pattern = f"{sample_id}*{path}.paths"
                matching_files = fnmatch.filter(os.listdir(paths_dir), pattern)

                # Print matching files
                if matching_files:
                    print(f"Matching files for {sample_id} with path {path}: {matching_files}")

                    for match in matching_files:
                        file_path = os.path.join(paths_dir, match)

                        # Open the matching file and read its contents
                        with open(file_path, 'r') as paths_file:
                            lines = paths_file.readlines()

                            # Check if the file has enough lines to append the 3rd and 4th lines
                            if len(lines) >= 4:
                                third_line = lines[2].strip()  # Get the 3rd line
                                fourth_line = lines[3].strip()  # Get the 4th line

                                # Append the 3rd and 4th lines to the current row
                                row.append(third_line)
                                row.append(fourth_line)

                            else:
                                print(f"File {file_path} does not have enough lines to append.")
                else:
                    print(f"No matching files found for {sample_id} with path {path}")

            # Write the modified or unmodified row to the new CSV
            csv_writer.writerow(row)

    print(f"Updated CSV saved as {updated_excel_doc}")

else:
    print(f"The file {excel_doc} does not exist.")
