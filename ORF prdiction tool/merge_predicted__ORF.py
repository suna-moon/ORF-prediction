import pandas as pd

day = [0, 4, 8]

for j in day:
    d0x = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{j}.txt', sep='\t')
    display(d0x)
    d0x = d0x[d0x['ORF_length']>=27] # stop codon을 포함하지 않으므로(기존 30에서 27로 변경)
    display(d0x)
    d0x.to_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{j}.length.txt', sep='\t', index=False)
