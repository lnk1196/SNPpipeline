# prep_seqs.py

import argparse
from subprocess import run

def run_snakemake(input_dir, sample_file=None, output_dir=None, cores=1):
    command = ["snakemake", "--use-conda", "--cores", str(cores)]

    if sample_file:
        command.extend(["--config", f"sample_file={sample_file}"])

    if output_dir:
        command.extend(["--config", f"output_dir={output_dir}"])

    run(command, cwd=input_dir, check=True)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run Snakemake workflow with specified inputs and optional outputs.")
    parser.add_argument("input_dir", help="Path to the input directory.")
    parser.add_argument("--sample_file", help="Path to the sample file (optional).", default=None)
    parser.add_argument("--output_dir", help="Path to the output directory (optional).", default=None)
    parser.add_argument("--cores", help="Number of cores to use for Snakemake (default is 1).", type=int, default=1)

    args = parser.parse_args()

    run_snakemake(args.input_dir, args.sample_file, args.output_dir, args.cores)
