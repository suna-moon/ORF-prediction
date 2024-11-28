# TE translation efficiency : RPF/RNA

# uORF TE : uORF RPF / uORF mRNA(CDS랑 동일)
# protein TE : CDS RPF / CDS mRNA
import pandas as pd


import matplotlib.pyplot as plt
import pandas as pd
import itertools

import numpy as np

from functools import reduce

# input = 'dORF'

# mRNA counts 중 필요한 것만 가져오기
mRNA = pd.read_csv('/Data_2/Jun/Adipocytes/counts/simple-features/RNA', sep='\t', skiprows=0, header=1)
# display(mRNA)
name_id = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/name_id.txt', sep='\t', header=0) # RNA counts에 필요한 gene id - gene name 정보
display(name_id)
mRNA_need = pd.merge(left=mRNA, right=name_id,
                how='right',
                left_on=['Geneid'],
                right_on=['gene_id']) # gene name - gene id merge
mRNA_need.columns = ['Geneid', 'Chr', 'Start', 'End', 'Strand', 'Length', 'gene_name', 'D0a', 'D0b', 'D0c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c', 'gene_id', 'gene_name2']

# mRNA RPKM 계산
df2 = mRNA_need.iloc[:, 7:16] # D0a ~ D8C
# display(df2)
# print(df2.sum(axis=0))
mRNAnormRPKMdf = df2.div(df2.sum(axis=0)) # total read count로 나누기
# display(mRNAnormRPKMdf)
mRNA_need = mRNA_need.set_index(mRNA_need['gene_name'])
mRNAnormRPKMdf = mRNAnormRPKMdf.set_index(mRNA_need['gene_name']) # Gene id

for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        mRNAnormRPKMdf[f'D{i}{j}'] = mRNAnormRPKMdf[f'D{i}{j}']/mRNA_need['Length']*1000000 # length 나누기
# print(mRNA_need['Length'])

# display(mRNAnormRPKMdf)

# log2
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        mRNAnormRPKMdf[f'D{i}{j}(log2)'] = np.log2(mRNAnormRPKMdf[f'D{i}{j}']) # log2(RPKM)
display(mRNAnormRPKMdf.iloc[:,9:18])


arr = [0, 4, 8]
for i, j in itertools.combinations(arr, 2):
    print(i, j)
    for k in ['a', 'b', 'c']:
        mRNAnormRPKMdf[f'D{j}{k}vsD{i}{k}'] = mRNAnormRPKMdf[f'D{j}{k}(log2)'] - mRNAnormRPKMdf[f'D{i}{k}(log2)'] # 비교, log2(a) - log2(b) 
display(mRNAnormRPKMdf.iloc[:, 18:27].sort_values('D4avsD0a', ascending=False))
# mRNAnormRPKMdf.to_csv(f'RPKM_log2_FC_RNA_{input}.txt', sep='\t')
# uORF TE (RPF/RNA)

uORF_rpkm = pd.read_csv('/Data_3/Suna/ORF/featureCounts/RPKM.txt', sep='\t', header=0, index_col='Geneid')

# dORF TE

# uORF_rpkm = pd.read_csv(f'/Data_3/Suna/ORF/featureCounts/RPKM_{input}.txt', sep='\t', header=0, index_col='Geneid')



# log(0) 처리
zero = []
for index, ser in uORF_rpkm.iterrows():
    for subindex, val in ser.items():
        if val == 0:
            zero.append(index)
            break
print(zero)
for i in zero:
    uORF_rpkm.loc[f'{i}'] = uORF_rpkm.loc[f'{i}'] + 1

# log2 씌우기
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        uORF_rpkm[f'D{i}{j}(log2)'] = np.log2(uORF_rpkm[f'D{i}{j}'])

# uORF_rpkm.to_csv(f'RPKM_uORF_RPF_log2_{input}.txt', sep='\t')

# FC
arr = [0, 4, 8]
for i, j in itertools.combinations(arr, 2):
    print(i, j)
    for k in ['a', 'b', 'c']:
        uORF_rpkm[f'D{j}{k}vsD{i}{k}'] = uORF_rpkm[f'D{j}{k}(log2)'] - uORF_rpkm[f'D{i}{k}(log2)'] # 비교, log2(a) - log2(b) 
display(uORF_rpkm.iloc[:, 18:27].sort_values('D4avsD0a', ascending=False))
# uORF_rpkm.to_csv(f'RPKM_uORF_log2_FC_RPF_{input}.txt', sep='\t')

# print(uORF_rpkm.iloc[:,18:27].columns) # 'D4avsD0a', 'D4bvsD0b', 'D4cvsD0c', 'D8avsD0a', 'D8bvsD0b', 'D8cvsD0c'

# TE = FC-FC
TE_uORF= pd.DataFrame()
for i, j in itertools.combinations(arr, 2):
    for k in ['a', 'b', 'c']:
        TE_uORF[f'D{j}{k}vsD{i}{k}'] = uORF_rpkm[f'D{j}{k}vsD{i}{k}'] - mRNAnormRPKMdf[f'D{j}{k}vsD{i}{k}'] 
display(TE_uORF)

# CDS TE
CDS_RPKM = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/RPKM/RPKM_CDS.txt', sep='\t', header=0, index_col='Geneid')
# CDS_RPKM = pd.read_csv(f'/Data_3/Suna/ORF/featureCounts/CDS_RPF/RPKM/RPKM_CDS_{input}.txt', sep='\t', header=0, index_col='Geneid')

display(CDS_RPKM)

zero = []
for index, ser in CDS_RPKM.iterrows():
    for subindex, val in ser.items():
        if val == 0:
            zero.append(index)
            break
print(zero)
for i in zero:
    CDS_RPKM.loc[f'{i}'] = CDS_RPKM.loc[f'{i}'] + 1

# log2
for i in range(0,12, 4):
    for j in ['a', 'b', 'c']:
        CDS_RPKM[f'D{i}{j}(log2)'] = np.log2(CDS_RPKM[f'D{i}{j}'])
display(CDS_RPKM)
# CDS_RPKM.to_csv(f'RPKM_CDS_RPF_log2_{input}.txt', sep='\t')

# FC
arr = [0, 4, 8]
for i, j in itertools.combinations(arr, 2):
    print(i, j)
    for k in ['a', 'b', 'c']:
        CDS_RPKM[f'D{j}{k}vsD{i}{k}'] = CDS_RPKM[f'D{j}{k}(log2)'] - CDS_RPKM[f'D{i}{k}(log2)'] # 비교, log2(a) - log2(b) 
display(CDS_RPKM.iloc[:, 18:27].sort_values('D4avsD0a', ascending=False))
# CDS_RPKM.to_csv(f'RPKM_CDS_log2_FC_RPF_{input}.txt', sep='\t')

# print(uORF_rpkm.iloc[:,18:27].columns) # 'D4avsD0a', 'D4bvsD0b', 'D4cvsD0c', 'D8avsD0a', 'D8bvsD0b', 'D8cvsD0c'

# TE = FC-FC
TE_CDS= pd.DataFrame()
for i, j in itertools.combinations(arr, 2):
    for k in ['a', 'b', 'c']:
        print(f'D{j}{k}vsD{i}{k}')
        TE_CDS[f'D{j}{k}vsD{i}{k}'] = CDS_RPKM[f'D{j}{k}vsD{i}{k}'] - mRNAnormRPKMdf[f'D{j}{k}vsD{i}{k}'] 
display(TE_CDS)

# TE_CDS, TE_uORF => heat map
df = pd.merge(left=TE_uORF.iloc[:, 0:6], right=TE_CDS.iloc[:, 0:6],
              left_index=True, right_index=True,
              suffixes=['_uORF', '_CDS'])

df['mean'] = df.iloc[:,3:6].mean(axis=1) # 정렬렬
df = df.sort_values('mean', ascending=False).iloc[:,0:12]
display(df)

# data
data_list = [df, uORF_rpkm.iloc[:, 18:24], CDS_RPKM.iloc[:, 18:24], mRNAnormRPKMdf.iloc[:, 18:24]]
df2 = pd.concat([df, uORF_rpkm.iloc[:, 18:24], CDS_RPKM.iloc[:, 18:24], mRNAnormRPKMdf.iloc[:, 18:24]], axis=1)
display(df2)
df2.to_csv('TE_uORF_data.txt', sep='\t')

# 정렬
# df2['mean'] = df2.iloc[:,0:3].mean(axis=1) # TE uORF D4/D0
df2['mean'] = df2.iloc[:,3:6].mean(axis=1) # TE uORF D8/D0
# df2['mean'] = df2.iloc[:,6:9].mean(axis=1) # TE CDS D4/D0
# df2['mean'] = df2.iloc[:,9:12].mean(axis=1) # TE CDS D8/D0
df2 = df2.sort_values('mean', ascending=False).iloc[:,0:30]

# df2['mean'] = df2.iloc[:,15:18].mean(axis=1) # RPKM FC uORF D8/D0
# df2 = df2.sort_values('mean', ascending=False).iloc[:,0:30]

# 그래프
# fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(12,8))

fig = plt.figure(figsize=(15,12)) # uORF
# fig = plt.figure(figsize=(15,3)) # dORF

gs = gridspec.GridSpec(nrows=1,
                       ncols=3,
                       width_ratios=[6, 6, 3]
                      )
ax1 = plt.subplot(gs[0])
# ax2 = plt.subplot(gs[1])
# ax3 = plt.subplot(gs[2])

# uORF
im = ax1.imshow(df2.iloc[:, 0:12], cmap='RdBu_r', aspect='auto', vmin=-3, vmax=3)
# im2 = ax2.imshow(df2.iloc[:,12:24], cmap='RdBu_r', aspect='auto', vmin=-8, vmax=8)
# im3 = ax3.imshow(df2.iloc[:, 24:30], cmap='RdBu_r', aspect='auto', vmin=-8, vmax=8)

# # dORF
# im = ax1.imshow(df2.iloc[:, 0:12], cmap='RdBu_r', aspect='auto', vmin=-3, vmax=3)
# im2 = ax2.imshow(df2.iloc[:,12:24], cmap='RdBu_r', aspect='auto', vmin=-5, vmax=5)
# im3 = ax3.imshow(df2.iloc[:, 24:30], cmap='RdBu_r', aspect='auto', vmin=-3, vmax=3)

# clolor bar
sh = 0.2 # uORF
# sh = 0.8 # dORF
cbar = ax1.figure.colorbar(im, ax=ax1, shrink=sh) 
cbar.set_label('TE ratio', size=10)
cbar.ax.tick_params(labelsize=7)
# cbar = ax2.figure.colorbar(im2, ax=ax2, shrink=sh) 
# cbar.set_label('RPKM FC', size=10)
# cbar.ax.tick_params(labelsize=7)
# cbar = ax3.figure.colorbar(im3, ax=ax3, shrink=sh) 
# cbar.set_label('RPKM FC', size=10)
# cbar.ax.tick_params(labelsize=7)

# 축
xlabel = ['D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c'] # 축
ax1.set_xticks(range(len(df2.iloc[:,12:24].columns)), labels=xlabel, fontsize=6)
ax1.set_yticks(range(len(df2.iloc[:,12:24].index)), labels=df2.iloc[:,12:24].index, fontsize=6)
ax1.tick_params(axis='x', length=0)
ax1.xaxis.tick_top() # x 축 위로

# ax2.set_xticks(range(len(df2.iloc[:,12:24].columns)), labels=xlabel, fontsize=6)
# ax2.set_yticks(range(len(df2.iloc[:,12:24].index)), labels=df2.iloc[:,12:24].index, fontsize=6)
# ax2.tick_params(axis='x', length=0)
# ax2.xaxis.tick_top()

# ax3.set_xticks(range(len(df2.iloc[:, 24:30].columns)), labels=['D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c'], fontsize=6)
# ax3.set_yticks(range(len(df2.iloc[:, 24:30].index)), labels=df2.iloc[:, 24:30].index, fontsize=6)
# ax3.tick_params(axis='x', length=0)
# ax3.xaxis.tick_top()

# 구분선
ax1.axvline(x=2.5, color='black', linewidth=0.7)
ax1.axvline(x=5.5, color='black', linewidth=1)
ax1.axvline(x=8.5, color='black', linewidth=0.7)

# ax2.axvline(x=2.5, color='black', linewidth=0.7)
# ax2.axvline(x=5.5, color='black', linewidth=1)
# ax2.axvline(x=8.5, color='black', linewidth=0.7)

# ax3.axvline(x=2.5, color='black', linewidth=1)

# 제목
ax1.set_title('TE \nRPF   CDS', fontsize=10)
# ax2.set_title('RPF', fontsize=10)
# ax3.set_title('RNA', fontsize=10)

# plt.title('uORF RPF', fontsize=10)

plt.tight_layout()
plt.savefig(f"TE_uORF_2.pdf", dpi = 300, bbox_inches = "tight")
plt.show()
