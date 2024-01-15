import os
from Bio import SeqIO

nucl_seqs_dir = '/mnt/scratch/brownlab/lkirsch/Ploidy/Felicity/nucl_matrix_constructor_out_Dec.08.2023/nucl_seqs'
gene_files = [os.path.join(nucl_seqs_dir, file) for file in os.listdir(nucl_seqs_dir) if file.endswith('.fas')]

org_dict = {}

for file in gene_files:
    gene = file.split('/')[-1].split('.fas')[0]
    with open(file, 'r') as infile:
        for record in SeqIO.parse(infile, 'fasta'):
            taxon = record.id  # Assuming the ID is the organism name
            # Create a new sequence name in the format organism_name_gene_name
            new_description = f"{taxon}_{gene}"
            record.id = new_description
            record.description = ""
            if taxon in org_dict:
                org_dict[taxon].append(record)
            else:
                org_dict[taxon] = [record]

# Create the "compiled_nucl" directory if it doesn't exist
output_dir = 'compiled_nucl'
os.makedirs(output_dir, exist_ok=True)

for k in org_dict:
    output_file_path = os.path.join(output_dir, f'{k}.fas')
    with open(output_file_path, 'w') as outfile:
        SeqIO.write(org_dict[k], outfile, 'fasta')
