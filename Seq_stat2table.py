##flagstatOut
#391043 + 0 in total (QC-passed reads + QC-failed reads)
#0 + 0 secondary
#935 + 0 supplementary
#0 + 0 duplicates
#282989 + 0 mapped (72.37% : N/A)
#390108 + 0 paired in sequencing
#195054 + 0 read1
#195054 + 0 read2
#263984 + 0 properly paired (67.67% : N/A)
#271364 + 0 with itself and mate mapped
#10690 + 0 singletons (2.74% : N/A)
#6300 + 0 with mate mapped to a different chr
#4872 + 0 with mate mapped to a different chr (mapQ>=5)

import os

file_list = os.listdir('.')
f_stat_list = [f_tmp for f_tmp in file_list if f_tmp.endswith('.flagstat.txt')]
f_stat_list_sort = sorted(f_stat_list)
f_out = open('200615_NextSeq_SH-SY5Y_RNASeq.SAMstat.txt','w')
f_out.write('Sample\tTotalReads\tMappedReads(%)\tPairedReads(%)\n')
for f_stat_nm in f_stat_list_sort:
    f_prefix = f_stat_nm.split('.')[0]
    f_stat = open(f_stat_nm,'r')
    num = 0
    for line in f_stat:
        tokens = line.strip().split(' ')
        num += 1
        if num == 1:
            total_read = tokens[0]
        elif num == 5:
            mapped_read = tokens[0]
            mapped_ratio = tokens[4].lstrip('(')
        elif num == 9:
            paired_read = tokens[0]
            paired_ratio = tokens[5].lstrip('(')

    f_out.write('%s\t%s\t%s(%s)\t%s(%s)\n'%(f_prefix,total_read,mapped_read,mapped_ratio,paired_read,paired_ratio))

f_out.close()
