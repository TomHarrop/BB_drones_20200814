#!/usr/bin/env python3

import csv
from Bio import SeqIO

# use the ragtag agp to orient the original scaffolds

agp_file = snakemake.input['agp']
fa_file = snakemake.input['fa']
out_file = snakemake.output['fa']

# dev
# agp_file = 'output/025_ragtag/BB31/ragtag.scaffolds.agp'
# fa_file = 'output/020_flye/BB31/assembly.fasta'
# out_file = 'test/ordered_contigs.fa'

with open(agp_file, 'rt') as f:
    my_csv = csv.reader(f, delimiter='\t')
    agp_lines = [x for x in my_csv if 
                 (not x[0].startswith("#") and
                  not x[4] == "U")]


fa_idx = SeqIO.index(fa_file, 'fasta')

outseqs = []

for line in agp_lines:
    print(line)
    record = fa_idx[line[5]]
    if line[8] == "-":
        record.seq = record.seq.reverse_complement()
    outseqs.append(record)

SeqIO.write(outseqs, out_file, "fasta")
