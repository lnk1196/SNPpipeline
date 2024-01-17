import sys
import glob
import os
import re

def oneline(db):
 infile = open(db, "r")
 lines = infile.readlines()
 infile.close()
 outfile = open(db,'w')
 for i,line in enumerate(lines):
    if line[0] == ('>'):
        if i>0:
            outfile.write("\n")
        outfile.write(line)
    else:
        #line = line.strip()
        #line = line.replace('-','')  ######## Add if you want to remove gaps!!!
        outfile.write(line.strip())
 outfile.close()




bacteria = ['atum', 'aaeo', 'aful', 'bant', 'bsui', 'bmal', 'bpse', 'cmaq', 'cjej', 'ckor', 'cpne', 'ctep', 'cbot', 'cper', 'cbur', 'deth', 'drad', 'ecol', 'ftul', 'halo', 'hwal', 'hbut', 'ihos', 'lmon', 'mgri', 'msed', 'msmi', 'mjan', 'mmar', 'mlep', 'mtub', 'nequ', 'nmar', 'rsol', 'rbal', 'rpro', 'rtyp', 'sent', 'sfle', 'saur', 'smar', 'spne', 'ssol', 'syne', 'tvol', 'tmar', 'tpal', 'vcho', 'wend', 'wsuc', 'ypes']

taxon = sys.argv[1]
infile = open("%s" %(taxon),"r")

oneline(sys.argv[2])

seqfile = open(sys.argv[2],"r") ### transcriptome fasta file NUCLEOTIDES
seqlines = seqfile.readlines()
seqfile.close()

'''
for line in seqlines:
    if line[0]==">":
        print (line)
'''

lines = infile.readlines()
infile.close()
outfile = open("%s-NoBact.out.fas" %(taxon),"w")

bactseqlist = []
for line in lines:
    #print (line)
    seq = line.split("\t")[0]
    #seq = seq.split("spades55_")[1]
    seq = seq.split(".")[0]
    group = line.split("\t")[1]
    taxonhit = line.split("\t")[2]
    taxonhit = taxonhit.split("|")[0]
    if taxonhit in bacteria:
        bactseqlist.append(seq)
        #print (taxonhit)
    else:
        pass

print ("DONE with FirstPass")
#print seqlines
#print (bactseqlist)
#print (len(bactseqlist))

for j in range(len(seqlines)):
    line =seqlines[j]
    #print line
    if line[0] == ">":
        header = line.split(">")[1]
        header = header.strip()
        header = header.split()[0] 
        #print (header)   
        #header = header.split('_i')[0]
        #print (header)
        if header in bactseqlist:
            pass
            print (header)
            print ("bacterial")
            #outfile.write(line)
            #next_l = lines[j+1]
            #outfile.write(next_l)
        else:
            #print header
            outfile.write(line)
            next_l = seqlines[j+1]
            outfile.write(next_l)
outfile.close()
