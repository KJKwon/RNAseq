import sys

file_num = int(sys.argv[1])
file_name_list = []
for num in range(file_num):
    tmp_f_name = input('File name?: ')
    file_name_list.append(tmp_f_name)

file_prefix = input('prefix?: ')

f_count_out = open(file_prefix+'_count.txt','w')
f_rpkm_out = open(file_prefix+'_rpkm.txt','w')
f_cpm_out = open(file_prefix+'_cpm.txt','w')

gene2count = dict()
gene2cpm = dict()
gene2rpkm = dict()
for f_name in file_name_list:
    f_input = open(f_name,'r')
    f_input.readline()
    f_input.readline()
    for line in f_input:
        token = line.strip().split('\t')
        gene_name = token[0]
        gene_count = token[1]
        gene_cpm = token[2]
        gene_rpkm = token[3]

        if gene_name not in gene2count:
            gene2count[gene_name] = []
        if gene_name not in gene2cpm:
            gene2cpm[gene_name] = []
        if gene_name not in gene2rpkm:
            gene2rpkm[gene_name] = []
        gene2count[gene_name].append(gene_count)
        gene2cpm[gene_name].append(gene_cpm)
        gene2rpkm[gene_name].append(gene_rpkm)
    
    f_input.close()

f_count_out.write('geneID\t')
f_cpm_out.write('geneID\t')
f_rpkm_out.write('geneID\t')
for num in range(file_num):
    f_count_out.write('%s_%s_count\t'%(file_prefix,num+1))
    f_cpm_out.write('%s_%s_cpm\t'%(file_prefix,num+1))
    f_rpkm_out.write('%s_%s_rpkm\t'%(file_prefix,num+1))

f_count_out.write('\n')
f_cpm_out.write('\n')
f_rpkm_out.write('\n')
for tmp_gene in sorted(gene2count.keys()):
    f_count_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2count[tmp_gene])))
    f_cpm_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2cpm[tmp_gene])))
    f_rpkm_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2rpkm[tmp_gene])))

f_count_out.close()
f_cpm_out.close()
f_rpkm_out.close()
