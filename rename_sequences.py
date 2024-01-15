import os
import argparse
import subprocess

def rename_seqs(input_dir, output_dir):
    # Check if the input directory exists
    if not os.path.exists(input_dir):
        print(f"Error: Directory {input_dir} does not exist.")
        exit(1)

    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read sample names from the sample file
    samples_file = os.path.join(input_dir, "samples.txt")
    with open(samples_file, 'r') as samples_file:
        samples = samples_file.read().splitlines()

    # Iterate through each file in the input directory
    for file_name in os.listdir(input_dir):
        file_path = os.path.join(input_dir, file_name)

        # Check if the file is a regular file
        if os.path.isfile(file_path):
            print(f"Processing file: {file_path}")

            # Extract the filename without the path
            filename = os.path.basename(file_path)

            # Create a temporary file for the processed content
            temp_file_path = os.path.join(output_dir, f"temp_{filename}")

            # Use awk to filter lines containing sample names and the following lines
            awk_command = f'awk -v samples="{",".join(samples)}" \'BEGIN {{split(samples, arr, ","); for (i in arr) sample_names[arr[i]] = 1;}} {{found = 0; for (sample in sample_names) {{if (index($0, sample) > 0) found = 1; break;}} if (found) {{sub(/_[0-9]+$/, "", $1); print $1; getline; print;}}}}\' "{file_path}" | awk -F\'_\' \'/^>/ {{print $1}} !/^>/ {{print}}\' > "{temp_file_path}"'
            subprocess.run(awk_command, shell=True, check=True)

            # Move the temporary file to the output directory
            processed_file_path = os.path.join(output_dir, filename)
            os.rename(temp_file_path, processed_file_path)

            print(f"Done processing file: {file_path}")

    print(f"Script completed successfully. Processed and modified files are in the {output_dir} directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename sequences in files based on sample names.")
    parser.add_argument("input_dir", help="Path to the input directory.")
    parser.add_argument("output_dir", help="Path to the output directory.")

    args = parser.parse_args()

    rename_seqs(args.input_dir, args.output_dir)
