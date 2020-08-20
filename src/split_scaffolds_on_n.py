#!/usr/bin/env python3

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import re

fasta_file = snakemake.input['fa']
map_file = snakemake.output['contig_map']
split_fasta = snakemake.output['contigs']

# dev
# fasta_file = 'output/025_ragtag/BB31/ragtag.scaffolds.fasta'
# map_file = 'test/split_contigs.csv'
# split_fasta = 'test/split_contigs.fa'

# initialise the map file
with open(map_file, 'wt') as f:
    f.write('old_id,new_id\n')

# initialise contig number
i = 0
outseqs = []

fasta = SeqIO.parse(fasta_file, 'fasta')
for myseq in fasta:
    my_splitseqs = re.compile("[nN]{100,}").split(str(myseq.seq))    
    for splitseq in my_splitseqs:
        i += 1
        my_newname = f'contig_{i:05}'
        # save the renaming details
        with open(map_file, 'a') as f:
            f.write(f'{myseq.id},{my_newname}\n')
        # generate a new Seq object
        outseqs.append(
            SeqRecord(
                Seq(splitseq),
                id=my_newname,
                name='',
                description=''))

SeqIO.write(outseqs, split_fasta, 'fasta')
