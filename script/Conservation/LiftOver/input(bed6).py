# liftover input file bed6 format : BED6: A BED file where each feature is described by chrom, start, end, name, score, and strand.
#                                                                           For example: chr1 11873 14409 uc001aaa.3 0 +
#                                                                                                                socre ? gray 
# input : cate(ORF_pos).csv,
# chrom : Col 10: strand_x
# start :Col 54: ORFstart => 0 based (-1 해야함)
# end : Col 55: ORFstop
# name : Col 7: gene_name
# score : exon 개수
# strand : Col 9: chrom_x


import pandas as pd

cate = ['uORF', 'dORF', 'lncRNA']


for c in cate:
    # position = pd.read_csv(f'/Data_3/Suna/ORF/mergeORF(day)/{i}(exon).csv', sep='\t')
    with open(f'/Data_3/Suna/ORF/mergeORF(day)/{c}(ORF_pos).csv', 'r') as f, open(f'/Data_3/Suna/ORF/Liftover/{c}.bed6', 'w') as out:
        write = []
        for line in f:
            line = line.rstrip('\n').split('\t')
            if line[0] == 'ORF_ID':
                continue
            start = line[53].split('.')
            stop = line[54].split('.')
            for l in range(len(start)):
                if line[9] == '-':
                    startW = int(stop[l]) -1
                    stopW = start[l]
                else:
                    startW = int(start[l]) -1  ### start 부분 -1 : 0 base 위함, 근데 -에서도 0 based가 start 쪽인가?의문
                    stopW = stop[l]
                write.extend([line[8], startW, stopW, line[6], 0, line[9]])
                out.write('\t'.join(map(str, write)) + '\n')
                write = []

# 중복 없애기

import pandas as pd

cate = ['uORF', 'dORF', 'lncRNA']

for i in cate:
    data = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/{i}.bed6', sep='\t', header=None)
    display(data)
    data = data.drop_duplicates(subset=[0, 1, 2], keep='first')
    display(data)
    data.to_csv(f'/Data_3/Suna/ORF/Liftover/{i}.bed6', sep='\t', header=False, index=False)


