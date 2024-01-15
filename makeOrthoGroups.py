import sys

# Define the list of bacterial taxons
bacteria = [
    'atum', 'aaeo', 'aful', 'bant', 'bsui', 'bmal', 'bpse', 'cmaq', 'cjej', 'ckor',
    'cpne', 'ctep', 'cbot', 'cper', 'cbur', 'deth', 'drad', 'ecol', 'ftul', 'halo',
    'hwal', 'hbut', 'ihos', 'lmon', 'mgri', 'msed', 'msmi', 'mjan', 'mmar', 'mlep',
    'mtub', 'nequ', 'nmar', 'rsol', 'rbal', 'rpro', 'rtyp', 'sent', 'sfle', 'saur',
    'smar', 'spne', 'ssol', 'syne', 'tvol', 'tmar', 'tpal', 'vcho', 'wend', 'wsuc',
    'ypes'
]

# Check for the correct number of command-line arguments
if len(sys.argv) != 3:
    print("Usage: python script.py <taxon> <genome>")
    sys.exit(1)

taxon = sys.argv[1]
genome = sys.argv[2]  # 'y' for genome, 'n' for transcriptome

# Read input file
with open(f"{taxon}/orthologGroups", "r") as infile:
    lines = infile.readlines()

# Initialize the list to store non-bacterial groups
group_list = []

# Process each line in the input file
for line in lines:
    parts = line.split("\t")
    group = parts[1]
    taxonhit = parts[2].split("|")[0]
    
    # Check if taxonhit is not in the list of bacteria
    if taxonhit not in bacteria:
        group_list.append(group)
    else:
        print(f"{group} hit bacteria")

# Write the output to a file
with open(f"{taxon}-NoBact.out", "w") as outfile:
    outfile.write(taxon + ",")

    for i in range(126536, 251276):
        group = f"OG5_{i}"
        if group in group_list:
            outfile.write("1")
        else:
            if genome == "y":
                outfile.write("0")
            else:
                outfile.write("-")
        
        outfile.write(",")

    outfile.write("\n")
