import os

TAXON = config['taxon']
READ1 = config['read1']
READ2 = config['read2']
THREADS = config['threads']

# Create a directory for the taxon
taxon_dir = os.path.join("output", TAXON)
os.makedirs(taxon_dir, exist_ok=True)

     
rule all:
    input:
        f'{TAXON}_summary.tsv'

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
        "yaml/bowtie2.yaml"
    log:
        f"{taxon_dir}/{TAXON}_ref_build.log"
    shell:
        '''
        echo "Organism Name: {TAXON}"
        bowtie2-build {input} {params.basename} > {taxon_dir}/{TAXON}_ref_build.log
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
    output: f"{taxon_dir}/{TAXON}_R.alignment.sam"
    threads: THREADS
    conda:
        "yaml/bowtie2.yaml"
    shell:
        '''
        bowtie2 -x {taxon_dir}/{TAXON}_restricted_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads}
        '''

rule SAM_to_BAM:
    input: f"{taxon_dir}/{TAXON}_R.alignment.sam"
    output: f"{taxon_dir}/{TAXON}_R.alignment.bam"
    conda:
        "yaml/samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''

rule SAM_sort:
    input: f"{taxon_dir}/{TAXON}_R.alignment.bam"
    output: f"{taxon_dir}/{TAXON}_R.alignment.sorted.bam"
    conda:
        "yaml/samtools.yaml"
    shell:
        '''
        samtools sort -o {output} {input}
        '''

rule mpileup:
    input: 
        f"restricted_assembly/{TAXON}.fas",
        f"{taxon_dir}/{TAXON}_R.alignment.sorted.bam"
    output: f"{taxon_dir}/{TAXON}_R.bcf"
    threads: THREADS 
    conda:
        "yaml/bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input: f"{taxon_dir}/{TAXON}_R.bcf"
    output: f"{taxon_dir}/{TAXON}_R.vcf"
    threads: THREADS 
    conda:
        "yaml/bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input:  f"{taxon_dir}/{TAXON}.vcf"
    output:  f"{TAXON}.snpden"
    conda:
        "yaml/vcftools.yaml"
    log:
        f"{TAXON}_snpden.log"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {TAXON} > {TAXON}_snpden.log 2>&1 
        '''
rule num_snps:
    input: f"{taxon_dir}/{TAXON}_ref_build.log"
    output: f"{taxon_dir}/{TAXON}_bases.txt" , f"{taxon_dir}/{TAXON}_taxon.txt"
    shell: 
        '''
        grep 'len:' {input} | head -n 1 | sed 's/len: //g' > {output[0]} 
        echo {TAXON} > {output[1]}
        '''

rule num_bases:
    input: f"{TAXON}_snpden.log"
    output: f"{taxon_dir}/{TAXON}_snps.txt"
    shell:
        '''
        grep -oP 'After filtering, kept \\K\\d+' {input} | sed -n '2p' > {output}
        '''

rule total_contigs:
    input: f"{taxon_dir}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    output: f"{taxon_dir}/{TAXON}_total_contigs.txt"
    shell: 
        '''
        grep -c ">" {input} > {output}
        '''

rule calc_SNPden:
    input: f'{TAXON}.snpden'
    output: f'{taxon_dir}/{TAXON}_calc.txt'
    shell:
        '''
        awk '{{ sum += $3 }} END {{ print sum / NR }}' {input} > {output}
        '''

rule unique_contigs:
    input: f'{TAXON}.snpden'
    output: f'{taxon_dir}/{TAXON}_IDs.txt'
    shell:
        '''
        awk '{{print $1}}' {input} | sort -u | wc -l > {output}
        '''

rule compile_data:
    input:
        taxon=f"{taxon_dir}/{TAXON}_taxon.txt",
        snps=f"{taxon_dir}/{TAXON}_snps.txt",
        bases=f"{taxon_dir}/{TAXON}_bases.txt",
        unique_contigs=f'{taxon_dir}/{TAXON}_IDs.txt',
        total_contigs=f"{taxon_dir}/{TAXON}_total_contigs.txt",
        calc_snpden=f'{taxon_dir}/{TAXON}_calc.txt'
    output: f'{TAXON}_summary.tsv'
    shell:
        '''
        paste {input} >> {output}
        '''