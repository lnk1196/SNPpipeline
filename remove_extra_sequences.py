import os
import argparse
import subprocess

def process_files(input_dir, output_dir=None):
    # Check if the input directory exists
    if not os.path.exists(input_dir):
        print(f"Error: Directory {input_dir} does not exist.")
        exit(1)

    # Create the output directory if it doesn't exist
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)

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

            # Use awk to remove lines containing "_q2c", "_q3c", "_q4c", or "_q5c" and the following lines,
            # and also skip lines with ".." in sequence names
            awk_command = f'awk \'/_q2c|_q3c|_q4c|_q5c/ {{flag=1; next}} flag {{flag=0; next}} /^>/ {{if ($0 ~ /\\.\\./) {{skip=1}} else {{skip=0}}}} !skip {{print}}\' "{file_path}" > "{temp_file_path}"'
            subprocess.run(awk_command, shell=True, check=True)

            # Move the temporary file to the output directory
            processed_file_path = os.path.join(output_dir, filename)
            os.rename(temp_file_path, processed_file_path)

            print(f"Done processing file: {file_path}")

    print(f"Script completed successfully. Processed files are in the {output_dir} directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process files in the input directory and save processed files in the output directory.")
    parser.add_argument("input_dir", help="Path to the input directory.")
    parser.add_argument("--output_dir", help="Path to the output directory (optional). If not provided, processed files will be saved in the input directory.")

    args = parser.parse_args()

    process_files(args.input_dir, args.output_dir)
