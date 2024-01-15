import os

TAXON = config['taxon']
READ1 = config['read1']
READ2 = config['read2']
ADAPTERS = config['adapters']

print("Organism Name: " + TAXON)

# Software paths
orthomcl_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/aa_seqs_OrthoMCL-5.fasta.removedspaces.fasta.dmnd"
bbmap_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/repair.sh"

# Python script paths
make_ortho_groups_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/makeOrthoGroups.py"
bacteria_check_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/bacCheck.py"
remove_bact_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/removeBac.py"

rule all:
    input:
        f"{TAXON}.snpden"

# Define rules
rule predict_orfs:
    input: f"raw_reads/{TAXON}.fas"
    output: f"{TAXON}.fas.transdecoder_dir/longest_orfs.pep",
            f"{TAXON}.fas.transdecoder_dir/longest_orfs.cds",
            f"{TAXON}.fas.transdecoder_dir/longest_orfs.gff3",
            f"{TAXON}.fas.transdecoder_dir/base_freqs.dat"
    threads: 16
    shell:
        '''
        echo "Organism Name: {TAXON}"
        /mnt/scratch/brownlab/lkirsch/Ploidy/scripts/TransDecoder.LongOrfs -t {input} -O f"{TAXON}.fas.transdecoder_dir"
        '''

rule rename_longest_orfs:
    input:  f"{TAXON}.fas.transdecoder_dir/longest_orfs.pep"
    output:  f"{TAXON}.faa"
    shell: "cp {input} {output}"

rule diamond_blast:
    input:  f"{TAXON}.faa"
    output:  f"{TAXON}.faa.OG5DMND.out"
    threads: 16
    conda:
        "diamond.yaml"
    shell:
        '''
        diamond blastp -f 6 -q {input} -p {threads} -o {output} -d {orthomcl_path} -k 1 -s sensitive
        '''

rule make_ortho_groups:
    input:  f"{TAXON}.faa.OG5DMND.out"
    output:  f"{TAXON}.faa.OG5DMND/orthologGroups"
    shell:
        '''
        python3 {make_ortho_groups_path} {input} 
        '''
        
rule bac_check:
    input:  f"{TAXON}.faa.OG5DMND"
    output:  f"{TAXON}.faa.OG5DMND-NoBact.out"
    shell:
        '''
        python3 {bacteria_check_path} {input} n
        '''
        
rule remove_bac:
    input:  f"{TAXON}.faa.OG5DMND/orthologGroups", f"raw_reads/{TAXON}.fas"
    output:  f"{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    shell:
        '''
        python3 {remove_bact_path} {input[0]} {input[1]} 
        '''

rule trimmomatic:
    input: 
        f"raw_reads/{READ1}", f"raw_reads/{READ2}"
    output:
         f'trimmed_{TAXON}_1.unpaired.fastq.gz',
         f"trimmed_{TAXON}_2.unpaired.fastq.gz",
         f'trimmed_{TAXON}_1.paired.fastq.gz',
         f"trimmed_{TAXON}_2.paired.fastq.gz"
    params:
        adapter = ADAPTERS
    threads: 16
    conda:
        "trimmomatic.yaml"
    shell:
        '''
        trimmomatic PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[2]} {output[1]} {output[3]} ILLUMINACLIP:/mnt/scratch/brownlab/lkirsch/Ploidy/adapters/{params.adapter}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
        '''

rule trimmed_1:
    input:  f"trimmed_{TAXON}_1.paired.fastq.gz",  f"trimmed_{TAXON}_1.unpaired.fastq.gz"
    output:  f"trimmed_{TAXON}_1.fastq.gz"
    shell:
        '''
        cat {input[0]} {input[1]} > {output}
        '''

rule trimmed_2:
    input:  f"trimmed_{TAXON}_2.paired.fastq.gz",  f"trimmed_{TAXON}_2.unpaired.fastq.gz"
    output:  f"trimmed_{TAXON}_2.fastq.gz"
    shell:
        '''
        cat {input[0]} {input[1]} > {output}
        '''

rule bbmap_repair:
    input:  f"trimmed_{TAXON}_1.fastq.gz",  f"trimmed_{TAXON}_2.fastq.gz"
    output:  f"fixed1_{TAXON}.fq.gz",  f"fixed2_{TAXON}.fq.gz",  f"singletons_{TAXON}.fq.gz"
    shell:
        '''
        {bbmap_path} in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} outsingle={output[2]}
        '''

rule bowtie_ref:
    input:  f"{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    params: 
        basename = f"{TAXON}_ref.fasta"
    output: 
         f"{TAXON}_ref.fasta.1.bt2",
         f"{TAXON}_ref.fasta.2.bt2",
         f"{TAXON}_ref.fasta.3.bt2",
         f"{TAXON}_ref.fasta.4.bt2",
         f"{TAXON}_ref.fasta.rev.1.bt2",
         f"{TAXON}_ref.fasta.rev.2.bt2"
    conda:
        "bowtie2.yaml"
    log:
        "{params.basename}.log"
    shell:
        '''
        bowtie2-build {input} {params.basename} > {params.basename}.log 2>&1
        '''

rule bowtie_align:
    input:
         f"{TAXON}_ref.fasta.1.bt2",
         f"{TAXON}_ref.fasta.2.bt2",
         f"{TAXON}_ref.fasta.3.bt2",
         f"{TAXON}_ref.fasta.4.bt2",
         f"{TAXON}_ref.fasta.rev.1.bt2",
         f"{TAXON}_ref.fasta.rev.2.bt2",
         f"fixed1_{TAXON}.fq.gz", 
         f"fixed2_{TAXON}.fq.gz"
    output:
         f"{TAXON}.alignment.sam"
    threads: 16
    log:
        f"{TAXON}_alignment.log"
    conda:
        "bowtie2.yaml"
    shell:
        '''
        bowtie2 -x {TAXON}_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads} > {TAXON}_alignment.log 2>&1
        '''

rule SAM_to_BAM:
    input:  f"{TAXON}.alignment.sam"
    output:  f"{TAXON}.alignment.bam"
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''

rule SAM_sort:
    input:  f"{TAXON}.alignment.bam"
    output:  f"{TAXON}.alignment.sorted.bam"
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools sort -o {output} {input}
        '''

rule mpileup:
    input: 
         f"{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas",
         f"{TAXON}.alignment.sorted.bam"
    output:  f"{TAXON}.bcf"
    threads: 16
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input:  f"{TAXON}.bcf"
    output:  f"{TAXON}.vcf"
    threads: 16
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input:  f"{TAXON}.vcf"
    output:  f"{TAXON}.snpden"
    conda:
        "vcftools.yaml"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {TAXON} > {output}.log
        '''

