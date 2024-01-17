import os

TAXON = config['taxon']
READ1 = config['read1']
READ2 = config['read2']
ADAPTERS = config['adapters']
#THREADS = config['threads']

print("Organism Name: " + TAXON)
#print(THREADS)

# Software paths
orthomcl_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/aa_seqs_OrthoMCL-5.fasta.removedspaces.fasta.dmnd"

# Python script paths
make_ortho_groups_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/makeOrthoGroups.py"
bacteria_check_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/bacCheck.py"
remove_bact_path = "/mnt/scratch/brownlab/lkirsch/Ploidy/scripts/removeBac.py"

rule all:
    input:
        f"{TAXON}.snpden"

# Define rules
rule trinity_assembly:
    input:
        f"raw_reads/{READ1}", f"raw_reads/{READ2}"
    output:
        f"{TAXON}_trinity_out_dir.Trinity.fasta",
        f"{TAXON}_trinity_out_dir/{READ1}.PwU.qtrim.fq",
        f"{TAXON}_trinity_out_dir/{READ2}.PwU.qtrim.fq"
    threads: 18 
    conda:
        "trinity.yaml"
    params:
        adapter = ADAPTERS
    shell:
        '''
        echo "Organism Name: {TAXON}"
        Trinity --seqType fq --max_memory 450G --left {input[0]} --right {input[1]} --CPU {threads} --trimmomatic --quality_trimming_params "ILLUMINACLIP:/mnt/scratch/brownlab/lkirsch/Ploidy/adapters/{params.ADAPTER}:2:30:10 SLIDINGWINDOW:4:5 LEADING:5 TRAILING:5 MINLEN:25" --output {output[0]}
        '''

rule rename_assembly:
    input: f"{TAXON}_trinity_out_dir.Trinity.fasta"
    output: f"{TAXON}.fas"
    shell:
        '''
        cp {input} {output}
        '''                                                                                                                                                                       

rule predict_orfs:
    input: f"{TAXON}.fas"
    output: f"{TAXON}.fas.transdecoder_dir/longest_orfs.pep",
            f"{TAXON}.fas.transdecoder_dir/longest_orfs.cds",
            f"{TAXON}.fas.transdecoder_dir/longest_orfs.gff3",
            f"{TAXON}.fas.transdecoder_dir/base_freqs.dat"
    threads: 18 
    conda:
        "transdecoder.yaml"
    shell:
        '''
        TransDecoder.LongOrfs -t {input}
        '''


rule rename_longest_orfs:
    input: f"{TAXON}.fas.transdecoder_dir/longest_orfs.pep"
    output: f"{TAXON}.faa"
    shell: "cp {input} {output}"

rule diamond_blast:
    input: f"{TAXON}.faa"
    output: f"{TAXON}.faa.OG5DMND.out"
    conda:
        "diamond.yaml"
    threads: 18 
    shell:
        '''
        diamond blastp -f 6 -q {input} -p {threads} -o {output} -d {orthomcl_path} -k 1 -s sensitive
        '''

rule make_ortho_groups:
    input: f"{TAXON}.faa.OG5DMND.out"
    output: f"{TAXON}.faa.OG5DMND/orthologGroups"
    shell:
        '''
        python3 {make_ortho_groups_path} {input} 
        '''
        
rule bac_check:
    input: f"{TAXON}.faa.OG5DMND"
    output: f"{TAXON}.faa.OG5DMND-NoBact.out"
    shell:
        '''
        python3 {bacteria_check_path} {input} n
        '''
        
rule remove_bac:
    input: f"{TAXON}.faa.OG5DMND/orthologGroups", f"{TAXON}.fas"
    output: f"{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    shell:
        '''
        python3 {remove_bact_path} {input[0]} {input[1]} 
        '''

rule bbmap_repair:
    input: f"{TAXON}_trinity_out_dir/{READ1}.PwU.qtrim.fq" , f"{TAXON}_trinity_out_dir/{READ2}.PwU.qtrim.fq" 
    output: f"fixed_1_{TAXON}.fq.gz", f"fixed_2_{TAXON}.fq.gz", f"singletons_{TAXON}.fq.gz"
    conda:
        "bbmap.yaml"
    shell:
        '''
        repair.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} outsingle={output[2]}
        '''

rule bowtie_ref:
    input: f"{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
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
    shell:
        '''
        bowtie2-build {input} {params.basename}
        '''

rule bowtie_align:
    input:
        f"{TAXON}_ref.fasta.1.bt2",
        f"{TAXON}_ref.fasta.2.bt2",
        f"{TAXON}_ref.fasta.3.bt2",
        f"{TAXON}_ref.fasta.4.bt2",
        f"{TAXON}_ref.fasta.rev.1.bt2",
        f"{TAXON}_ref.fasta.rev.2.bt2",
        f"fixed_1_{TAXON}.fq.gz", 
        f"fixed_2_{TAXON}.fq.gz"
    output:
        f"{TAXON}.alignment.sam"
    conda:
        "bowtie2.yaml"
    threads: 18 
    shell:
        '''
        bowtie2 -x {TAXON}_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads}
        '''


rule SAM_to_BAM:
    input: f"{TAXON}.alignment.sam"
    output: f"{TAXON}.alignment.bam"
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''


rule SAM_sort:
    input: f"{TAXON}.alignment.bam"
    output: f"{TAXON}.alignment.sorted.bam"
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
    output: f"{TAXON}.bcf"
    threads: 18 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input: f"{TAXON}.bcf"
    output: f"{TAXON}.vcf"
    threads: 18 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input: f"{TAXON}.vcf"
    output: f"{TAXON}.snpden"
    conda:
        "vcftools.yaml"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {TAXON} > {output}.log
        '''