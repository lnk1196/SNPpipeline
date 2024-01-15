import argparse
import os

parser = argparse.ArgumentParser(description="Run Snakemake pipeline with Trinity assemble or not")
parser.add_argument("samples", help="path to samples.txt")
args = parser.parse_args()

with open(args.samples) as f:
    for line in f:
        taxon, read1, read2, adapters, assembly = [x.strip() for x in line.split()]
        if assembly == "Y":
            os.system(f"snakemake -s SNP_w_assembly.smk --cores 16 --use-conda --config taxon={taxon} read1={read1} read2={read2} adapters={adapters}")
        elif assembly == "N":
            os.system(f"snakemake -s SNP_wo_assembly.smk --cores 16 --use-conda --config taxon={taxon} read1={read1} read2={read2} adapters={adapters}")
        else:
            print(f"Invalid value for assembly in line: {line}")
