import argparse
import os
import shutil

# Parse the command-line arguments
parser = argparse.ArgumentParser(description="Copy and rename longest_orf.cds files.")
parser.add_argument("sample_file", help="Path to the sample file (toRUN2.txt)")
parser.add_argument("target_dir", help="Destination directory for copied files")
args = parser.parse_args()

# Create the target directory if it doesn't exist
os.makedirs(args.target_dir, exist_ok=True)

# Read the taxon names from the sample file
with open(args.sample_file, "r") as file:
    taxon_names = [line.split()[0] for line in file.readlines()]

# Process each taxon
for taxon in taxon_names:
    source_file = f"{taxon}.fas.transdecoder_dir/longest_orfs.pep"
    target_file = f"{args.target_dir}/{taxon}.faa"

    # Check if the source file exists
    if os.path.exists(source_file):
        try:
            # Copy and rename the file
            shutil.copy(source_file, target_file)
            print(f"Copied and renamed {source_file} to {target_file}")
        except OSError as error:
            print(f"Error copying file: {error}")
    else:
        print(f"Source file not found: {source_file}")
