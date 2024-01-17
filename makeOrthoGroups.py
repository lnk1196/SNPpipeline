import os
import sys
import glob

'''
#!/bin/bash
#MOAB/Torque submission script for Shadow
#PBS -l nodes=1:ppn=12,walltime=48:00:00
#PBS -N Mal249
cd $PBS_O_WORKDIR
blastp -num_threads 12 -max_target_seqs 1 -query Mal249.faa -outfmt "6 std qcovhsp" -db /scratch/mbrown/Backup/work/mbrown/DB/OrthoMCL/aa_seqs_OrthoMCL-5.fasta.removedspaces.fasta -out Mal249.out 
'''

input=open(sys.argv[1], 'r')
lines=input.readlines()
input.close()
outfile = open("orthologGroups","w")
#outfile = open(sys.argv[2],"w")

dir = sys.argv[1].split(".out")[0]

try:
 os.system('rm -r %s' %(dir))
except:
 pass

os.system('mkdir %s' %(dir))

def lines_seen(input):
 intemp = open(input, "r")
 templines = intemp.readlines()
 intemp.close()
 tempout = open(input, "w")
 lines_seen = set() 
 
 for line in templines:
    geneid = line.split()[0]
    if geneid not in lines_seen:
        tempout.write(line)
        lines_seen.add(geneid)
 tempout.close() 


for line in lines:
    #qcovhsp = float(line.split("\t")[12])
    qseqid = line.split("\t")[0]
    evalue = float(line.split("\t")[10])
    sseqid = line.split("\t")[1]
    try:    
        oid = line.split("OG5_")[1]
        oid = oid.split("|")[0]
        #if qcovhsp >= 50:
        if evalue <= 0.00001:
                outfile.write(qseqid)
                outfile.write("\t")
                outfile.write("OG5_")
                outfile.write(oid)
                outfile.write("\t")
                outfile.write(sseqid)
                outfile.write("\t")
                outfile.write(str(evalue))
                outfile.write("\n")
        #   else:
        #       pass
        else:
            pass
    except:
        pass

outfile.close()

lines_seen("orthologGroups")

os.system('mv orthologGroups %s' %(dir))


