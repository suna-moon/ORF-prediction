tri = ['a', 'b', 'c']
day = [0, 4, 8]

for i in day:
    for j in tri:


        with open(f'/Data_3/Suna/RiboCode/output/D{i}{j}/RiboCode_ORFs_result_collapsed.txt', 'r') as f, open(f'D{i}{j}.txt', 'w') as out:
            for read in f:
                final =[]
                if read.split('\t')[29] ==  'adjusted_pval':
                    tapsplit = read.rstrip('\n').split('\t')
                    out.write('\t'.join(map(str, tapsplit)) + '\n')
                elif float(read.split('\t')[29]) < 0.05:
                    tapsplit = read.rstrip('\n').split('\t')
                    out.write('\t'.join(map(str, tapsplit)) + '\n')
