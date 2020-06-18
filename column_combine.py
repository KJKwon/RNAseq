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

title_list = []
gene2count = dict()
gene2cpm = dict()
gene2rpkm = dict()
for f_name in file_name_list:
    f_input = open(f_name,'r')
    f_input.readline()
    f_input.readline()
    title_list.append(f_name.split('.count+cpm+rpkm.txt')[0])
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

f_count_out.write('geneID\t%s\n'%('\t'.join(title_list)))
f_cpm_out.write('geneID\t%s\n'%('\t'.join(title_list)))
f_rpkm_out.write('geneID\t%s\n'%('\t'.join(title_list)))

gene_list = sorted(gene2count.keys())
for tmp_gene in gene_list:
    f_count_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2count[tmp_gene])))
    f_cpm_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2cpm[tmp_gene])))
    f_rpkm_out.write('%s\t%s\n'%(tmp_gene,'\t'.join(gene2rpkm[tmp_gene])))

f_count_out.close()
f_cpm_out.close()
f_rpkm_out.close()
