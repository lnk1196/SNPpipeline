import os

TAXON = config['taxon']
READ1 = config['read1']
READ2 = config['read2']
#THREADS = config['threads']

# Create a directory for the taxon
taxon_dir = os.path.join("output", TAXON)
os.makedirs(taxon_dir, exist_ok=True)

     
rule all:
    input:
        expand(os.path.join(taxon_dir, f"{TAXON}_R.snpden"))

# define rules
rule bowtie_ref:
    input: f"restricted_assembly/{TAXON}.fas"
    output:
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.1.bt2",
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.2.bt2",
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.3.bt2",
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.4.bt2",
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.rev.1.bt2",
        f"{taxon_dir}/{TAXON}_restricted_ref.fasta.rev.2.bt2"
    params: 
        basename = f"{taxon_dir}/{TAXON}_restricted_ref.fasta"
    conda:
        "bowtie2.yaml"
    shell:
        '''
        echo "Organism Name: {TAXON}"
        bowtie2-build {input} {params.basename}
        '''


rule bowtie_align:
    input:
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.1.bt2",
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.2.bt2",
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.3.bt2",
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.4.bt2",
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.rev.1.bt2",
            f"{taxon_dir}/{TAXON}_restricted_ref.fasta.rev.2.bt2",
            f"fixed_rawreads/{READ1}", 
            f"fixed_rawreads/{READ2}"
    output: os.path.join(taxon_dir, f"{TAXON}_R.alignment.sam")
    threads: 65
    conda:
        "bowtie2.yaml"
    shell:
        '''
        bowtie2 -x {taxon_dir}/{TAXON}_restricted_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads}
        '''

rule SAM_to_BAM:
    input: os.path.join(taxon_dir, f"{TAXON}_R.alignment.sam")
    output: os.path.join(taxon_dir, f"{TAXON}_R.alignment.bam")
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''

rule SAM_sort:
    input: os.path.join(taxon_dir, f"{TAXON}_R.alignment.bam")
    output: os.path.join(taxon_dir, f"{TAXON}_R.alignment.sorted.bam")
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools sort -o {output} {input}
        '''

rule mpileup:
    input: 
            os.path.join("restricted_assembly", f"{TAXON}.fas"),
            os.path.join(taxon_dir, f"{TAXON}_R.alignment.sorted.bam")
    output: os.path.join(taxon_dir, f"{TAXON}_R.bcf")
    threads: 65 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input: os.path.join(taxon_dir, f"{TAXON}_R.bcf")
    output: os.path.join(taxon_dir, f"{TAXON}_R.vcf")
    threads: 65 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input: os.path.join(taxon_dir, f"{TAXON}_R.vcf")
    output: os.path.join(taxon_dir, f"{TAXON}_R.snpden")
    conda:
        "vcftools.yaml"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {taxon_dir}/{TAXON}_R > {output}.log
        '''
