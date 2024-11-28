import pandas as pd
import seaborn as sns
import numpy as np

input = 'uORF'

df = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/counts_output/featureCountesCDS_output.txt', sep='\t', skiprows=1)
# df = pd.read_csv(f'/Data_3/Suna/ORF/featureCounts/CDS_RPF/counts_output/featureCountesCDS_output_{input}.txt', sep='\t', skiprows=1)

df.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'D0a', 'D0b', 'D0c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c']
display(df)
df2 = df.iloc[:, 6:15]

# display(df2)
# display(df2.sum(axis=0))
normRPKMdf = df2.div(df2.sum(axis=0)) # total read count로 나누기
# display(normRPKMdf)
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        normRPKMdf[f'D{i}{j}'] = normRPKMdf[f'D{i}{j}']/df['Length']*1000000 # length 나누기
normRPKMdf = normRPKMdf.set_index(df['Geneid']) # Gene id

# display(df['Length'])

display(normRPKMdf)
normRPKMdf.to_csv(f'RPKM_CDS.txt', sep='\t')
# normRPKMdf.to_csv(f'RPKM_CDS_{input}.txt', sep='\t')
