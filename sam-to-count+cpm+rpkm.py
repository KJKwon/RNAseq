import sys
import os
import re

f_sam_best = sys.argv[1]
f_sam_best_prefix = f_sam_best.split('.')[0]

sam_best = open(f_sam_best,'r')
gene2len = dict()
gene2pair_count = dict()

for line in sam_best:
    if line.startswith('@SQ'):
       token = line.strip().split()
       gene_ID = re.sub('^SN:','',token[1])
       gene_length = re.sub('^LN:','',token[2])
       gene2len[gene_ID] = int(gene_length)
       continue

    token = line.strip().split('\t')
    if len(token) < 6:
        continue

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
        if match_gene_ID not in gene2pair_count:
            gene2pair_count[match_gene_ID] = 0
        gene2pair_count[match_gene_ID] += 1

sam_best.close()

total_pair_count = sum(gene2pair_count.values())
f_output = open('%s.count+cpm+rpkm.txt'%(f_sam_best_prefix),'w')
f_output.write('#PairCount : %d\n'%(total_pair_count))
f_output.write('GeneID\tPairCount\tPairCPM\tPairRPKM\tGeneLength\n')
for t_ID in sorted(gene2len.keys()):
    t_len = gene2len[t_ID]
    if t_ID in gene2pair_count:
        t_pair_count = float(gene2pair_count[t_ID])
        t_CPM = t_pair_count/(total_pair_count/1000000.00)
        t_len = gene2len[t_ID]
        t_RPKM = t_pair_count/((total_pair_count*t_len)/1000000000.00)
        f_output.write('%s\t%d\t%.2f\t%.2f\t%d\n'%(t_ID,int(t_pair_count),t_CPM,t_RPKM,t_len))
    else:
        f_output.write('%s\t0\t0\t0\t%d\n'%(t_ID, t_len))

f_output.close()
