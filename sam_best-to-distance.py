import sys

##sam_file
#FS10000436:28:BPG80716-1926:1:1101:2330:1000    83      Chn1|ENSRNOT00000068745 2634    48      87S29M35S       =       2631    -32     CCCCCCCCCCCCTTTTTTCAAGCAGAAGACGGCATACGAGATTCCTCTACGTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGTATACACACATACACACCACACACATACCACACACACTAATTCCACAAACCTCAATAAATTCCTATA FFF,FFFFFF:,FF::FF:FFF,FFFFFF:FFFFF:F:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF: NM:i:1  MD:Z:22C6       MC:Z:32M119S    AS:i:24 XS:i:20 XA:Z:AABR07049768.1|ENSRNOT00000088532,-845,101S20M30S,0;AABR07049768.1|ENSRNOT00000088532,-463,101S20M30S,0;
f_in_nm = sys.argv[1]
f_in = open(f_in_nm,'r')
f_out_nm = f_in_nm.split('.')[0]+'_paired_distance.txt'
f_out = open(f_out_nm,'w')
for line in f_in:
    if line.startswith('@'):
        continue
    else:
        token = line.strip().split('\t')
        header = token[0]
        starts = int(token[3])
        ends = int(token[7])
        if starts > ends:
            continue
        else:
            read_length = len(token[9])
            distance = (ends+read_length) - starts + 1
            if distance < 0:
                distance = 0
            f_out.write('%s\t%d\n'%(header, distance))

f_out.close()
