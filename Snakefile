#!/usr/bin/env python3

import snakemake


###########
# GLOBALS #
###########

porechop = 'shub://TomHarrop/ont-containers:porechop_0.2.4'
bbmap = 'shub://TomHarrop/seq-utils:bbmap_38.76'


def combine_indiv_reads(wildcards):
    my_read_path = f'data/{wildcards.indiv}/pass/{{read}}.fastq.gz'
    my_output_path = f'output/010_porechop/{wildcards.indiv}.fastq'
    my_read_names = snakemake.io.glob_wildcards(my_read_names).read
    my_output = snakemake.io.expand(my_output_path, read=my_read_names)
    return(sorted(set(my_output)))


#########
# RULES #
#########

rule target:
    input:
        expand('output/010_porechop/{indiv}.fastq',
               indiv=['BB31', 'BB55'])

rule combine_indiv_reads:
    input:
        combine_indiv_reads
    output:
        'output/010_porechop/{indiv}.fastq'
    singularity:
        bbmap
    shell:
        'cat {input} > {output}'

rule porechop:
    input:
        'data/{indiv}/pass/{read}.fastq.gz'
    output:
        'output/010_porechop/{indiv}/{read}.fastq'
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

