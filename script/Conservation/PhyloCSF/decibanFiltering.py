
import pandas as pd


input = 'lncRNA' # dORF, uORF

df = pd.read_csv(f'/Data_3/Suna/ORF/PhyloCSF/{input}.txt', sep='\t', header=None)
display(df)
filename = []

filterdf = df[df[2] > 0]
display(filterdf)
for i in filterdf[0].values:
    filename.append(i.split('/')[8].split('_')[0])
filterdf[0] = filename
filterdf.to_csv(f'{input}filter0.txt', header=None, index=None, sep='\t')
