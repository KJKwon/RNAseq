import sys
import os
import re

f_sam = sys.argv[1]
f_sam_prefix = f_sam.split('.')[0]

sam_file = open(f_sam,'r')
gene2len = dict()
gene2contig = dict()
total_read = 0
best_read = 0
f_sam_best_nm = f_sam.split('.sam')[0] +'.best.sam'
f_sam_best = open(f_sam_best_nm,'w')
for line in sam_file:
    if line.startswith('@SQ'):
        f_sam_best.write(line)
        continue
    token = line.strip().split('\t')

    if len(token) < 6:
        continue

    total_read += 1
    match_gene_ID = token[2]
    flag = int(token[1])
    mapq = int(token[4])
    seq1_tag = token[2]
    seq2_tag = token[6]
    seq1_pos = int(token[3])
    seq2_pos = int(token[7])

    if(seq1_tag == '*' or flag & 4 ):
        continue
    if(mapq <= 0):
        continue
    if seq1_pos == seq2_pos:
        continue
    if seq1_tag == seq2_tag or seq2_tag == '=':
        best_read += 1
        f_sam_best.write(line)

f_sam_best.close()
sam_file.close()
