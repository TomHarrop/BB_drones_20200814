#!/usr/bin/env python3

import csv
import logging
from Bio import SeqIO

# use the ragtag agp to orient the original scaffolds

agp_file = snakemake.input['agp']
fa_file = snakemake.input['fa']
out_file = snakemake.output['fa']
log_file = snakemake.log[0]

# dev
# agp_file = 'output/025_ragtag/BB31/ragtag.scaffolds.agp'
# fa_file = 'output/020_flye/BB31/assembly.fasta'
# out_file = 'test/ordered_contigs.fa'
# log_file = 'test/orient.log'


logging.basicConfig(
    filename=log_file,
    format='%(asctime)s %(levelname)-8s %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO)

logging.info(f'Parsing {agp_file}')
with open(agp_file, 'rt') as f:
    my_csv = csv.reader(f, delimiter='\t')
    agp_lines = [x for x in my_csv if
                 (not x[0].startswith('#') and
                  not x[4] == 'U')]

logging.info(f'Indexing {fa_file}')
fa_idx = SeqIO.index(fa_file, 'fasta')

outseqs = []

logging.info(f'Generating output sequences')
for line in agp_lines:
    logging.info(line)
    record = fa_idx[line[5]]
    if line[8] == '-':
        record.seq = record.seq.reverse_complement()
    outseqs.append(record)

logging.info(f'Writing output sequences to {out_file}')
SeqIO.write(outseqs, out_file, 'fasta')
