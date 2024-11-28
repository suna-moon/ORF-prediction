# %%
# matagene data 생성

import pysam
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

dataStart = []
dataStop = []
readLength = []
stoppickgenes = []
tri = ['a', 'b', 'c']
time = ['0', '4', '8']
# D0b.rep.bam

def count_elements(lst):
    counts = {}
    for item in lst:
        if item in counts:
            counts[item] += 1
        else:
            counts[item] = 1
    return counts

df = pd.read_csv('data.txt', sep='\t')
# filepath = '/Data_2/Jun/Adipocytes/tr-aln/novaseq/D0a.rep.bam'

tri = ['a', 'b', 'c']
time = ['0', '4', '8']

lenlist =df[['5UTR_length', 'CDS_length', '3UTR_length']].values.tolist()
geneName = df['transcript_id'].values.tolist()
datadict = dict(zip(geneName, lenlist))
for i in time:
    for j in tri:
        filepath = f'/Data_2/Jun/Adipocytes/tr-aln/novaseq/D{i}{j}.rep.bam'
        bamfile = pysam.AlignmentFile(filepath, 'rb')
        dataStart = []
        dataStop = []

        for line in tqdm(bamfile):
            total = line.to_string().split('\t')

            if line.get_tag('NH') >= 2:   # unique read
                continue
            if datadict.get(total[2], False):
                if datadict[total[2]][0] == 0 or datadict[total[2]][1] == 0:
                    continue
                else :
                    a = int(total[3]) - datadict[total[2]][0] - 1
                    b = int(total[3]) - ( datadict[total[2]][0] + datadict[total[2]][1] - 2 ) 
                    dataStart.append(a)
                    dataStop.append(b)
             
        resultdictStart = count_elements(dataStart)
        resultdictStop = count_elements(dataStop)

        dfStart = pd.DataFrame({key:[val] for key, val in resultdictStart.items()}).T.sort_index().loc[-100:100, :]
        dfStop = pd.DataFrame({key:[val] for key, val in resultdictStop.items()}).T.sort_index().loc[-100:100, :]
        # a = dfStart.sum(axis=0) + dfStop.sum(axis=0)

        dfStart.to_csv(f'/Data_3/Suna/start/D{i}{j}.csv')
        dfStop.to_csv(f'/Data_3/Suna/stop/D{i}{j}.csv')







