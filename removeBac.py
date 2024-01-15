import sys
import os

def remove_gaps_from_fasta(fasta_file):
    with open(fasta_file, "r") as infile, open(fasta_file + '.tmp', 'w') as outfile:
        in_sequence = False
        for line in infile:
            if line.startswith('>'):
                if in_sequence:
                    outfile.write('\n')
                outfile.write(line)
                in_sequence = True
            else:
                line = line.strip().replace('-', '')  # Remove gaps
                outfile.write(line)
        if in_sequence:
            outfile.write('\n')

def extract_bacterial_sequences(taxon_file, sequence_file, output_file, bacteria_list):
    with open(taxon_file, "r") as infile:
        taxon_lines = infile.readlines()

    with open(sequence_file, "r") as seqfile:
        sequence_lines = seqfile.readlines()

    bact_seq_set = set()
    for line in taxon_lines:
        parts = line.strip().split("\t")
        taxonhit = parts[2].split("|")[0]
        if taxonhit in bacteria_list:
            bact_seq_set.add(parts[0])

    with open(output_file, "w") as outfile:
        in_sequence = False
        for line in sequence_lines:
            if line.startswith('>'):
                header = line.split(">")[1].split()[0].strip()
                if header in bact_seq_set:
                    in_sequence = True
                else:
                    in_sequence = False
                    outfile.write(line)
            elif in_sequence:
                outfile.write(line)

if __name__ == '__main__':
    if len(sys.argv) != 4:
        print("Usage: python script.py <taxon_file> <sequence_file> <output_file>")
        sys.exit(1)

    taxon_file = sys.argv[1]
    sequence_file = sys.argv[2]
    output_file = sys.argv[3]

    bacteria = ['atum', 'aaeo', 'aful', 'bant', 'bsui', 'bmal', 'bpse', 'cmaq', 'cjej', 'ckor', 'cpne', 'ctep', 'cbot', 'cper', 'cbur', 'deth', 'drad', 'ecol', 'ftul', 'halo', 'hwal', 'hbut', 'ihos', 'lmon', 'mgri', 'msed', 'msmi', 'mjan', 'mmar', 'mlep', 'mtub', 'nequ', 'nmar', 'rsol', 'rbal', 'rpro', 'rtyp', 'sent', 'sfle', 'saur', 'smar', 'spne', 'ssol', 'syne', 'tvol', 'tmar', 'tpal', 'vcho', 'wend', 'wsuc', 'ypes']

    remove_gaps_from_fasta(sequence_file)
    extract_bacterial_sequences(taxon_file, sequence_file + '.tmp', output_file, bacteria)

    os.remove(sequence_file + '.tmp')
