#!/usr/bin/env python3

import snakemake
import tempfile


#############
# FUNCTIONS #
#############

def combine_indiv_reads(wildcards):
    my_read_path = f'data/{wildcards.indiv}/pass/{{read}}.fastq.gz'
    my_output_path = f'output/010_porechop/{wildcards.indiv}/{{read}}.fastq'
    my_read_names = snakemake.io.glob_wildcards(my_read_path).read
    my_output = snakemake.io.expand(my_output_path, read=my_read_names)
    return(sorted(set(my_output)))


###########
# GLOBALS #
###########

bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.86'
flye = 'shub://TomHarrop/assemblers:flye_2.8'
minimap = 'shub://TomHarrop/singularity-containers:minimap2_2.17r941'
porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
samtools = 'shub://TomHarrop/align-utils:samtools_1.10'
sniffles = 'shub://TomHarrop/variant-utils:sniffles_53b7500'

indivs = ['BB31', 'BB55']

#########
# RULES #
#########

wildcard_constraints:
    indiv = '|'.join(indivs)

rule target:
    input:
        expand('output/{folder}/{indiv}.sorted.bam.bai',
               indiv=indivs,
               folder=['040_wga']),
        expand('output/050_sniffles/{indiv}.vcf.gz',
               indiv=indivs)


# SVs
rule sniffles:
    input:
        'output/030_mapped/{indiv}.sorted.bam'
    output:
        'output/050_sniffles/{indiv}.vcf'
    log:
        'output/logs/sniffles.{indiv}.log'
    threads:
        min(workflow.cores, 64)
    singularity:
        sniffles
    shell:
        'sniffles '
        '-m {input} '
        '-v {output} '
        f'--tmp_file {tempfile.mkdtemp()} '
        '-t {threads} '
        '&> {log}'

# WGA
rule wga:
    input:
        fa = 'output/020_flye/{indiv}/assembly.fasta',
        ref = 'output/000_ref/ref.mmi'
    output:
        pipe('output/040_wga/{indiv}.sam')
    log:
        'output/logs/wga.{indiv}.log'
    threads:
        min(workflow.cores, 64)
    singularity:
        minimap
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax asm5 '
        '--MD '
        '{input.ref} '
        '{input.fa} '
        '> {output} '
        '2> {log}'


# REFERENCE MAP
rule sort:
    input:
        'output/{folder}/{indiv}.bam'
    output:
        'output/{folder}/{indiv}.sorted.bam'
    log:
        'output/logs/sort.{folder}.{indiv}.log'
    threads:
        min(workflow.cores, 4)
    singularity:
        samtools
    shell:
        'samtools sort '
        '-@ {threads} '
        '-l 9 '
        '{input} '
        '> {output} '
        '2> {log}'

rule sam_to_bam:
    input:
        'output/{folder}/{indiv}.sam'
    output:
        pipe('output/{folder}/{indiv}.bam')
    log:
        'output/logs/sam_to_bam.{folder}.{indiv}.log'
    threads:
        1
    singularity:
        samtools
    shell:
        'samtools view -bh -u {input} '
        '>> {output} '
        '2> {log}'

rule map_to_genome:
    input:
        fq = 'output/010_porechop/{indiv}.fastq.gz',
        ref = 'output/000_ref/ref.mmi'
    output:
        pipe('output/030_mapped/{indiv}.sam')
    params:
        rg = '\'@RG\\tID:{indiv}\\tSM:{indiv}\''
    log:
        'output/logs/map_to_genome.{indiv}.log'
    threads:
        min(workflow.cores, 64)
    singularity:
        minimap
    shell:
        'minimap2 '
        '-t {threads} '
        '-ax map-ont '
        '--MD '
        '-R {params.rg} '
        '{input.ref} '
        '{input.fq} '
        '> {output} '
        '2> {log}'

rule prepare_ref:
    input:
        'data/GCF_003254395.2_Amel_HAv3.1_genomic.fna'
    output:
        'output/000_ref/ref.mmi'
    log:
        'output/logs/prepare_ref.log'
    threads:
        3
    singularity:
        minimap
    shell:
        'minimap2 '
        '-x map-ont '
        '-d {output} '
        '{input} '
        '2> {log}'


# DE NOVO ASSEMBLY
rule flye:
    input:
        fq = 'output/010_porechop/{indiv}.fastq.gz'
    output:
        'output/020_flye/{indiv}/assembly.fasta'
    params:
        outdir = 'output/020_flye/{indiv}',
        size = '250m'
    threads:
        min(128, workflow.cores)
    log:
        'output/logs/flye.{indiv}.log'
    singularity:
        flye
    shell:
        'flye '
        # '--resume '
        '--nano-raw {input.fq} '
        '--genome-size {params.size} '
        '--out-dir {params.outdir} '
        '--threads {threads} '
        '&>> {log}'

# PROCESS READS
rule combine_indiv_reads:
    input:
        combine_indiv_reads
    output:
        pipe('output/010_porechop/{indiv}.fastq')
    singularity:
        bbmap
    shell:
        'cat {input} > {output}'

rule porechop:
    input:
        'data/{indiv}/pass/{read}.fastq.gz'
    output:
        temp('output/010_porechop/{indiv}/{read}.fastq')
    log:
        'output/logs/porechop.{indiv}.{read}.log'
    threads:
        1
    singularity:
        porechop
    shell:
        'porechop '
        '-i {input} '
        '-o {output} '
        '--verbosity 1 '
        '--threads {threads} '
        '--discard_middle '
        '&> {log}'

# GENERICS
rule compress_reads:
    input:
        'output/{path}/{file}.fastq'
    output:
        'output/{path}/{file}.fastq.gz'
    threads:
        min(workflow.cores, 10)
    singularity:
        bbmap
    shell:
        'pigz -c --best {input} > {output}'


rule index_bamfile:
    input:
        'output/{folder}/{indiv}.sorted.bam'
    output:
        'output/{folder}/{indiv}.sorted.bam.bai'
    log:
        'output/logs/index_bamfile.{folder}.{indiv}.log'
    singularity:
        samtools
    shell:
        'samtools index {input} 2> {log}'

rule index_vcf:
    input:
        'output/{folder}/{file}.vcf'
    output:
        gz = 'output/{folder}/{file}.vcf.gz',
        tbi = 'output/{folder}/{file}.vcf.gz.tbi'
    singularity:
        samtools
    shell:
        'bgzip -c {input} > {output.gz} '
        '; '
        'tabix -p vcf {output.gz}'
