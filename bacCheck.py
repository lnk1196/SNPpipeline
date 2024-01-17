import sys
import glob
import os

bacteria = ['atum', 'aaeo', 'aful', 'bant', 'bsui', 'bmal', 'bpse', 'cmaq', 'cjej', 'ckor', 'cpne', 'ctep', 'cbot', 'cper', 'cbur', 'deth', 'drad', 'ecol', 'ftul', 'halo', 'hwal', 'hbut', 'ihos', 'lmon', 'mgri', 'msed', 'msmi', 'mjan', 'mmar', 'mlep', 'mtub', 'nequ', 'nmar', 'rsol', 'rbal', 'rpro', 'rtyp', 'sent', 'sfle', 'saur', 'smar', 'spne', 'ssol', 'syne', 'tvol', 'tmar', 'tpal', 'vcho', 'wend', 'wsuc', 'ypes']

taxon = sys.argv[1]
genome = sys.argv[2] ### y or n == y for genome, n for transcriptome
infile = open("%s/orthologGroups" %(taxon),"r")
lines = infile.readlines()
infile.close()
outfile = open("%s-NoBact.out" %(taxon),"w")

list = []
for line in lines:
    group = line.split("\t")[1]
    taxonhit = line.split("\t")[2]
    taxonhit = taxonhit.split("|")[0]
    if taxonhit not in bacteria:
        list.append(group)
    else:
        pass
        print ("%s hit bacteria" %(group))

outfile.write(taxon)
outfile.write(",")

if genome == "n":
 for i in range(126536,251276):
    group = "OG5_%s" %(i)
    if group in list:
        outfile.write("1")
        outfile.write(",")
    else:
        outfile.write("-")
        outfile.write(",")

elif genome == "y":
 for i in range(126536,251276):
    group = "OG5_%s" %(i)
    if group in list:
        outfile.write("1")
        outfile.write(",")
    else:
        outfile.write("0")
        outfile.write(",")

outfile.write("\n")
outfile.close()

