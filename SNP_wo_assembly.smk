import os

TAXON = config['taxon']
READ1 = config['read1']
READ2 = config['read2']
ADAPTERS = config['adapters']

output_directory = f"output/{TAXON}"

print("Organism Name: " + TAXON)

# Software paths
orthomcl_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/aa_seqs_OrthoMCL-5.fasta.removedspaces.fasta.dmnd"

# Python script paths
make_ortho_groups_path = "makeOrthoGroups.py"
bacteria_check_path = "bacCheck.py"
remove_bact_path = "removeBac.py"

rule all:
    input:
        f"{TAXON}.snpden"

# Define rules
rule predict_orfs:
    input: f"raw_reads/{TAXON}.fas"
    output: directory(f"{output_directory}/{TAXON}.fas.transdecoder_dir"),
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.pep",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.cds",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.gff3",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/base_freqs.dat"
    threads: 40
    conda:
        "yaml/transdecoder.yaml"
    shell:
        '''
        echo "Organism Name: {TAXON}"
        TransDecoder.LongOrfs -t {input} -O "{output_directory}/{TAXON}.fas.transdecoder_dir"
        '''

rule rename_longest_orfs:
    input:  f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.pep"
    output:  f"{output_directory}/{TAXON}.faa"
    shell: "cp {input} {output}"

rule diamond_blast:
    input:  f"{output_directory}/{TAXON}.faa"
    output:  f"{output_directory}/{TAXON}.faa.OG5DMND.out"
    threads: 40
    conda:
        "yaml/diamond.yaml"
    shell:
        '''
        diamond blastp -f 6 -q {input} -p {threads} -o {output} -d {orthomcl_path} -k 1 -s sensitive
        '''

rule make_ortho_groups:
    input:  f"{output_directory}/{TAXON}.faa.OG5DMND.out"
    output:  f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups"
    shell:
        '''
        python3 {make_ortho_groups_path} {input} 
        '''
        
rule bac_check:
    input:  f"{output_directory}/{TAXON}.faa.OG5DMND"
    output:  f"{output_directory}/{TAXON}.faa.OG5DMND-NoBact.out"
    shell:
        '''
        python3 {bacteria_check_path} {input} n
        '''
        
rule remove_bac:
    input:  f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups", f"raw_reads/{TAXON}.fas"
    output:  f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    shell:
        '''
        python3 {remove_bact_path} {input[0]} {input[1]} 
        '''

rule trimmomatic:
    input: 
        f"raw_reads/{READ1}", f"raw_reads/{READ2}"
    output:
         f'{output_directory}/trimmed_{TAXON}_1.unpaired.fastq.gz',
         f"{output_directory}/trimmed_{TAXON}_2.unpaired.fastq.gz",
         f'{output_directory}/trimmed_{TAXON}_1.paired.fastq.gz',
         f"{output_directory}/trimmed_{TAXON}_2.paired.fastq.gz"
    params:
        adapter = ADAPTERS
    threads: 40
    conda:
        "yaml/trimmomatic.yaml"
    shell:
        '''
        trimmomatic PE -threads {threads} {input[0]} {input[1]} {output[0]} {output[2]} {output[1]} {output[3]} ILLUMINACLIP:/mnt/scratch/brownlab/lkirsch/Ploidy/adapters/{params.adapter}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25
        '''

rule trimmed_1:
    input:  f"{output_directory}/trimmed_{TAXON}_1.paired.fastq.gz",  f"{output_directory}/trimmed_{TAXON}_1.unpaired.fastq.gz"
    output:  f"{output_directory}/trimmed_{TAXON}_1.fastq.gz"
    shell:
        '''
        cat {input[0]} {input[1]} > {output}
        rm {output_directory}/trimmed_{TAXON}_1.paired.fastq.gz {output_directory}/trimmed_{TAXON}_1.unpaired.fastq.gz
        '''

rule trimmed_2:
    input:  f"{output_directory}/trimmed_{TAXON}_2.paired.fastq.gz",  f"{output_directory}/trimmed_{TAXON}_2.unpaired.fastq.gz"
    output:  f"{output_directory}/trimmed_{TAXON}_2.fastq.gz"
    shell:
        '''
        cat {input[0]} {input[1]} > {output}
        rm {output_directory}/trimmed_{TAXON}_2.paired.fastq.gz {output_directory}/trimmed_{TAXON}_2.unpaired.fastq.gz
        '''

rule bbmap_repair:
    input:  f"{output_directory}/trimmed_{TAXON}_1.fastq.gz",  f"{output_directory}/trimmed_{TAXON}_2.fastq.gz"
    output:  f"{output_directory}/fixed_1_{TAXON}.fq.gz",  f"{output_directory}/fixed_2_{TAXON}.fq.gz",  f"{output_directory}/singletons_{TAXON}.fq.gz"
    conda:
        "yaml/bbmap.yaml"
    shell:
        '''
        repair.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} outsingle={output[2]}
        '''

rule bowtie_ref:
    input:  f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    params: 
        basename = f"{output_directory}/{TAXON}_ref.fasta"
    output: 
         f"{output_directory}/{TAXON}_ref.fasta.1.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.2.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.3.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.4.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.rev.1.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.rev.2.bt2"
    conda:
        "yaml/bowtie2.yaml"
    log:
        f"{output_directory}/{TAXON}_ref_build.log"
    shell:
        '''
        bowtie2-build {input} {params.basename} > {output_directory}/{TAXON}_ref_build.log
        '''

rule bowtie_align:
    input:
         f"{output_directory}/{TAXON}_ref.fasta.1.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.2.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.3.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.4.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.rev.1.bt2",
         f"{output_directory}/{TAXON}_ref.fasta.rev.2.bt2",
         f"{output_directory}/fixed_1_{TAXON}.fq.gz", 
         f"{output_directory}/fixed_2_{TAXON}.fq.gz"
    output:
         f"{output_directory}/{TAXON}.alignment.sam"
    threads: 40
    log:
        f"{output_directory}/{TAXON}_alignment.log"
    conda:
        "yaml/bowtie2.yaml"
    shell:
        '''
        bowtie2 -x {output_directory}/{TAXON}_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads} > {output_directory}/{TAXON}_alignment.log
        '''

rule SAM_to_BAM:
    input:  f"{output_directory}/{TAXON}.alignment.sam"
    output:  f"{output_directory}/{TAXON}.alignment.bam"
    conda:
        "yaml/samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''

rule SAM_sort:
    input:  f"{output_directory}/{TAXON}.alignment.bam"
    output:  f"{output_directory}/{TAXON}.alignment.sorted.bam"
    conda:
        "yaml/samtools.yaml"
    shell:
        '''
        samtools sort -o {output} {input}
        rm {output_directory}/{TAXON}.alignment.sam
        '''

rule mpileup:
    input: 
         f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas",
         f"{output_directory}/{TAXON}.alignment.sorted.bam"
    output:  f"{output_directory}/{TAXON}.bcf"
    threads: 40
    conda:
        "yaml/bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input:  f"{output_directory}/{TAXON}.bcf"
    output:  f"{output_directory}/{TAXON}.vcf"
    threads: 40
    conda:
        "yaml/bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input:  f"{output_directory}/{TAXON}.vcf"
    output:  f"{TAXON}.snpden"
    conda:
        "yaml/vcftools.yaml"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {TAXON} > {output}.log
        '''

