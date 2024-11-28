# 중복 없애기

import pandas as pd
input = 'lncRNA' #uORF, dORF

# data = pd.read_csv(f'/Data_3/Suna/ORF/mergeORF(day)/uORF(ORF_pos).csv', sep='\t', header=None)
data = pd.read_csv(f'/Data_3/Suna/ORF/mergeORF(day)/{input}(ORF_pos).csv', sep='\t', header=None)

data = data.drop_duplicates(subset=[0, 1, 2], keep='first')
# data.to_csv('/Data_3/Suna/ORF/featureCounts/uORF(ORF_pos_norep).txt', sep='\t', header=False, index=False)
data.to_csv(f'/Data_3/Suna/ORF/featureCounts/{input}(ORF_pos_norep).txt', sep='\t', header=False, index=False)

# PhyloCSF score > 0인 uORF의 chromosome cordinate 불러오기
# phylocsf 결과 /Data_3/Suna/ORF/PhyloCSF/uORFfilter0.txt
# chromosome cordinate 정보 /Data_3/Suna/ORF/featureCounts/uORF(ORF_pos_norep).txt

# phyloDF = pd.read_csv('/Data_3/Suna/ORF/PhyloCSF/uORFfilter0.txt', sep='\t', header=None)
phyloDF = pd.read_csv(f'/Data_3/Suna/ORF/PhyloCSF/{input}filter0.txt', sep='\t', header=None)
# cordinateDF = pd.read_csv('/Data_3/Suna/ORF/featureCounts/uORF(ORF_pos_norep).txt', sep='\t')
cordinateDF = pd.read_csv(f'/Data_3/Suna/ORF/featureCounts/{input}(ORF_pos_norep).txt', sep='\t')

cordinateDF = cordinateDF[['gene_name', 'chrom_x', 'ORFstart', 'ORFstop', 'strand_x']]
mergeDF = pd.merge(phyloDF, cordinateDF, how='left', left_on=0, right_on='gene_name')
mergeDF.rename(columns={2:'score(decibans)'}, inplace=True)
# mergeDF.iloc[:,2:8].to_csv(f'uORFcordinate.txt', sep='\t', index=None)
mergeDF.iloc[:,2:8].to_csv(f'{input}cordinate.txt', sep='\t', index=None)
mergeDF = mergeDF.iloc[:,2:8]
display(mergeDF)

with open(f'/Data_3/Suna/ORF/featureCounts/{input}cordinate.txt', 'r') as f, open(f'/Data_3/Suna/ORF/featureCounts/{input}cordinate.saf', 'w') as out:
    for i in f:
        item = i.rstrip('\n').split('\t')
        if item[0] == 'score(decibans)':
            write = 'GeneID','chr', 'Start', 'End','Strand'
            out.write('\t'.join(write) + '\n')
            continue
        if len(item[3].split('.'))>=2:
            for j, k in zip(item[3].split('.'), item[4].split('.')):
                if item[5] == '+':
                    write = '\t'.join(item[1:3])+'\t'+ j+'\t'+k + '\t' + str(item[5])
                    # if j > k :
                    #     print(write)
                    out.write(write + '\n')
                else:
                    write = '\t'.join(item[1:3])+'\t'+ k+'\t'+j + '\t' + str(item[5])
                    # if k > j:
                    #     print(write)
                    out.write(write + '\n')

        elif item[5] == '-':
                write = '\t'.join(item[1:3])+'\t'+ item[4]+'\t'+item[3] + '\t' + str(item[5])
                # if item[4] > item[3]:
                #     print(write)
                out.write(write +'\n')
        else:
            write = item[1:6]
            out.write('\t'.join(write) + '\n')
