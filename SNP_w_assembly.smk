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
rule trinity_assembly:
    input:
        f"raw_reads/{READ1}", f"raw_reads/{READ2}"
    output:
        f"{output_directory}/{TAXON}_trinity_out_dir.Trinity.fasta",
        f"{output_directory}/{TAXON}_trinity_out_dir/{READ1}.PwU.qtrim.fq",
        f"{output_directory}/{TAXON}_trinity_out_dir/{READ2}.PwU.qtrim.fq"
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
    input: f"{output_directory}/{TAXON}_trinity_out_dir.Trinity.fasta"
    output: f"{output_directory}/{TAXON}.fas"
    shell:
        '''
        cp {input} {output}
        '''                                                                                                                                                                       

rule predict_orfs:
    input: f"{output_directory}/{TAXON}.fas"
    output: f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.pep",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.cds",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.gff3",
            f"{output_directory}/{TAXON}.fas.transdecoder_dir/base_freqs.dat"
    threads: 18 
    conda:
        "transdecoder.yaml"
    shell:
        '''
        TransDecoder.LongOrfs -t {input}
        '''


rule rename_longest_orfs:
    input: f"{output_directory}/{TAXON}.fas.transdecoder_dir/longest_orfs.pep"
    output: f"{output_directory}/{TAXON}.faa"
    shell: "cp {input} {output}"

rule diamond_blast:
    input: f"{output_directory}/{TAXON}.faa"
    output: f"{output_directory}/{TAXON}.faa.OG5DMND.out"
    conda:
        "diamond.yaml"
    threads: 18 
    shell:
        '''
        diamond blastp -f 6 -q {input} -p {threads} -o {output} -d {orthomcl_path} -k 1 -s sensitive
        '''

rule make_ortho_groups:
    input: f"{output_directory}/{TAXON}.faa.OG5DMND.out"
    output: f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups"
    shell:
        '''
        python3 {make_ortho_groups_path} {input} 
        '''
        
rule bac_check:
    input: f"{output_directory}/{TAXON}.faa.OG5DMND"
    output: f"{output_directory}/{TAXON}.faa.OG5DMND-NoBact.out"
    shell:
        '''
        python3 {bacteria_check_path} {input} n
        '''
        
rule remove_bac:
    input: f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups", f"{output_directory}/{TAXON}.fas"
    output: f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    shell:
        '''
        python3 {remove_bact_path} {input[0]} {input[1]} 
        '''

rule bbmap_repair:
    input: 
        f"{output_directory}/{TAXON}_trinity_out_dir/{READ1}.PwU.qtrim.fq" , 
        f"{output_directory}/{TAXON}_trinity_out_dir/{READ2}.PwU.qtrim.fq" 
    output: 
        f"{output_directory}/fixed_1_{TAXON}.fq.gz", 
        f"{output_directory}/fixed_2_{TAXON}.fq.gz", 
        f"{output_directory}/singletons_{TAXON}.fq.gz"
    conda:
        "bbmap.yaml"
    shell:
        '''
        repair.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} outsingle={output[2]}
        '''

rule bowtie_ref:
    input: f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
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
        "bowtie2.yaml"
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
    conda:
        "bowtie2.yaml"
    log:
        f"{output_directory}/{TAXON}.alignment.log"
    threads: 18 
    shell:
        '''
        bowtie2 -x {TAXON}_ref.fasta -1 {input[6]} -2 {input[7]} -S {output} -p {threads} > {output_directory}/{TAXON}_alignment.log
        '''


rule SAM_to_BAM:
    input: f"{output_directory}/{TAXON}.alignment.sam"
    output: f"{output_directory}/{TAXON}.alignment.bam"
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools view -S -b {input} > {output}
        '''


rule SAM_sort:
    input: f"{output_directory}/{TAXON}.alignment.bam"
    output: f"{output_directory}/{TAXON}.alignment.sorted.bam"
    conda:
        "samtools.yaml"
    shell:
        '''
        samtools sort -o {output} {input}
        rm {output_directory}/{TAXON}.alignment.sam
        '''


rule mpileup:
    input: 
        f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas",
        f"{output_directory}/{TAXON}.alignment.sorted.bam"
    output: f"{output_directory}/{TAXON}.bcf"
    threads: 18 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools mpileup -O b -o {output} -f {input[0]} --threads {threads} -q 20 -Q 20 {input[1]}
        '''

rule call_SNPs:
    input: f"{output_directory}/{TAXON}.bcf"
    output: f"{output_directory}/{TAXON}.vcf"
    threads: 18 
    conda:
        "bcftools.yaml"
    shell:
        '''
        bcftools call --skip-variants indels --multiallelic-caller --variants-only -O v {input} -o {output} --threads {threads}
        '''

rule SNP_density:
    input: f"{output_directory}/{TAXON}.vcf"
    output: f"{TAXON}.snpden"
    conda:
        "vcftools.yaml"
    shell:
        '''
        vcftools --vcf {input} --SNPdensity 1000 --out {TAXON} > {output}.log
        '''

rule num_snps:
    input: f"{output_directory}/{TAXON}_ref_build.log"
    output: f"{output_directory}/{TAXON}_snps.txt" , f"{output_directory}/{TAXON}_taxon.txt"
    shell: 
        '''
        grep 'len:' {input} | head -n 1 | sed 's/len: //g' > {output[0]} 
        echo {TAXON} > {output[1]}
        '''

rule num_bases:
    input: f"{TAXON}_snpden.log"
    output: f"{output_directory}/{TAXON}_bases.txt"
    shell:
        '''
        grep -oP 'After filtering, kept \K\d+' {input} | sed -n '2p' > {output}
        '''

rule total_contigs:
    input: f"{output_directory}/{TAXON}.faa.OG5DMND/orthologGroups-NoBact.out.fas"
    output: f"{output_directory}/{TAXON}_total_contigs.txt"
    shell: "grep -c ">" {input} > {output}"

rule calc_SNPden:
    input: f'{TAXON}.snpden'
    output: f'{output_directory}/{TAXON}_calc.txt'
    shell:
        '''
        awk '{ sum += $3 } END { print sum / NR }' {input} > {output}
        '''

rule unique_contigs:
    input: f'{TAXON}.snpden'
    output: f'{output_directory}/{TAXON}_IDs.txt'
    shell:
        '''
        awk '{print $1}' {input} | sort -u | wc -l >> {output}
        '''

rule compile_data:
    input:
        taxon=f"{output_directory}/{TAXON}_taxon.txt"
        snps=f"{output_directory}/{TAXON}_snps.txt",
        bases=f"{output_directory}/{TAXON}_bases.txt",
        unique_contigs=f'{output_directory}/{TAXON}_IDs.txt',
        total_contigs=f"{output_directory}/{TAXON}_total_contigs.txt",
        calc_snpden=f'{output_directory}/{TAXON}_calc.txt'
    output: f'{TAXON}_summary.tsv'
    shell:
        '''
        paste {input} >> {output}
        '''