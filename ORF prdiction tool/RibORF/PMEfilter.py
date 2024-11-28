# %%
tri = ['a', 'b', 'c']
day = [0, 4, 8]

orfidindex = 'TranscriptID:chromatin:strand|RankNumber|transcriptLength:startCodonPosition:stopCodonPosition|candidateORFType|startCodonType'


for i in day:
    for j in tri:


        with open(f'/Data_3/Suna/RibORFfile/result/D{i}{j}/repre.valid.pred.pvalue.parameters.txt', 'r') as f, open(f'D{i}{j}.txt', 'w') as out:
        # with open('/Data_3/Suna/RibORFfile/result/D0a/repre.valid.pred.pvalue.parameters.txt', 'r') as f:
            for read in f:
                final =[]
                if read.split('\t')[12] ==  'PME':
                    tapsplit = read.rstrip('\n').split('\t')
                    final.append(orfidindex.split('|')[0].split(':')[0])
                    for k in orfidindex.split('|')[1:4]:
                        final.append(k)
                    for l in tapsplit[1::]:
                        final.append(l)
                    # print([orfidindex.split('|')[0].split(':')[0], orfidindex.split('|')[1:4], tapsplit])
                    out.write('\t'.join(map(str, final)) + '\n')
                    # out.write(read)
                elif float(read.split('\t')[12]) >= 0.6:
                    tapsplit = read.rstrip('\n').split('\t')
                    final.append(tapsplit[0].split('|')[0].split(':')[0])
                    for k in tapsplit[0].split('|')[1:4]:
                        final.append(k)
                    for l in tapsplit[1::]:
                        final.append(l)
                    out.write('\t'.join(map(str, final)) + '\n')




