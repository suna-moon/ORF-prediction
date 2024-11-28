import pandas as pd

tri = ['a', 'b', 'c']
day = [0, 4, 8]

for j in day:
    d0 = pd.DataFrame()
    for i in tri:
        d0x = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/p-value_filter/D{j}{i}.txt', sep='\t')
        d0 = pd.concat([d0, d0x])
        if i == 'b' or i == 'c':
            d0 = d0[d0.duplicated(subset=['ORF_ID'])]
    d0.to_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{j}.txt', sep='\t', index=False)
