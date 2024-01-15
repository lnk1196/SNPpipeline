import os
import sys

def remove_duplicate_lines(input_file, output_file):
    lines_seen = set()
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for line in infile:
            gene_id = line.split()[0]
            if gene_id not in lines_seen:
                outfile.write(line)
                lines_seen.add(gene_id)

def process_blast_results(input_file, output_dir):
    with open(input_file, 'r') as infile, open('orthologGroups', 'w') as outfile:
        for line in infile:
            fields = line.split('\t')
            qseqid = fields[0]
            evalue = float(fields[10])
            sseqid = fields[1]
            try:
                oid = fields[1].split('OG5_')[1].split('|')[0]
                if evalue <= 0.00001:
                    outfile.write(f"{qseqid}\tOG5_{oid}\t{sseqid}\t{evalue}\n")
            except:
                pass
    remove_duplicate_lines('orthologGroups', f'{output_dir}/orthologGroups')

if __name__ == '__main__':
    if len(sys.argv) != 3:
        print('Usage: python script.py <input_file> <output_directory>')
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    process_blast_results(input_file, output_dir)
