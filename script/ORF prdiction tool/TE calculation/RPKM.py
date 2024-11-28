import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import math
import numpy as np
import itertools

input = 'lncRNA' # 'lncRNA' uORF
# df = pd.read_csv(f'/Data_3/Suna/ORF/featureCounts/featureCountes_output_{input}.txt', sep='\t', skiprows=1)
df = pd.read_csv('/Data_3/Suna/ORF/featureCounts/featureCountes_output.txt', sep='\t', skiprows=1)
df.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'D0a', 'D0b', 'D0c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c']
display(df)
df2 = df.iloc[:, 6:15]

print(df2.sum(axis=0))
normRPKMdf = df2.div(df2.sum(axis=0)) # total read count로 나누기
# normRPKMdf = df2.copy(deep=True) # lncRNA만

display(normRPKMdf)
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        normRPKMdf[f'D{i}{j}'] = normRPKMdf[f'D{i}{j}']/df['Length']*1000000 # length 나누기
print(df['Length'])
normRPKMdf = normRPKMdf.set_index(df['Geneid']) # Gene id
# display(df['Length'])

display(normRPKMdf)
normRPKMdf.to_csv(f'RPKM_{input}.txt', sep='\t')

normRPKMdf.loc['Hr'] = normRPKMdf.loc['Hr'] + 1
for i in range(0,12, 4):
    
    for j in ['a', 'b', 'c']:
        normRPKMdf[f'D{i}{j}(log2)'] = np.log2(normRPKMdf[f'D{i}{j}']) # log2(RPKM)
display(normRPKMdf.iloc[:,9:18])

# 비교 빼기
arr = [0, 4, 8]
for i, j in itertools.combinations(arr, 2):
    print(i, j)
    for k in ['a', 'b', 'c']:
        normRPKMdf[f'D{j}{k}vsD{i}{k}'] = normRPKMdf[f'D{j}{k}(log2)'] - normRPKMdf[f'D{i}{k}(log2)'] # 비교, log2(a) - log2(b) 
display(normRPKMdf.iloc[:, 18:27].sort_values('D4avsD0a', ascending=False))


df2 = df.set_index('Geneid').iloc[:, 5:15]
display(df2)
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        normRPKMdf[f'D{i}{j}'] = df2[f'D{i}{j}']

writeDF = df2.merge(normRPKMdf.iloc[:, 9:27], left_on='Geneid', right_on='Geneid')
display(normRPKMdf.iloc[:, 9:27])
display(writeDF)
writeDF.to_csv(f'RPKM_FC_{input}.txt', sep='\t')
