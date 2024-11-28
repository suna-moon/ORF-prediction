import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import math
from scipy import stats

# x축 data : /Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/TE_uORF_data.txt
# TE, RPF, RNA
# 0:12, 12:24, 24:36
TE_FC = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/TE_uORF_data.txt', sep='\t', header = 0)
# display(TE)
TE_FC = TE_FC.iloc[:, 0:7]
# display(TE)
TE_FC['meanD4'] = TE_FC.iloc[:, 1:4].mean(axis=1)
TE_FC['meanD8'] = TE_FC.iloc[:, 4:7].mean(axis=1)
TE_FC = TE_FC.set_index('Unnamed: 0')
display(TE_FC)

# y 축 가공전
# RPKM RNA /Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_log2_RNA.txt
rpkm_RPF = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_uORF_RPF_log2.txt', sep='\t', header=0, index_col='Geneid')
rpkm_RPF = rpkm_RPF.iloc[:, 0:9]
display(rpkm_RPF)
# RPKM RPF(uORF) 
rpkm_RNA = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_log2_RNA.txt', sep='\t', header = 0, index_col='gene_name')
rpkm_RNA = rpkm_RNA.iloc[:, 0:9]
display(rpkm_RNA)

TE = rpkm_RPF/rpkm_RNA
display(TE)

# paired t-test
p = pd.DataFrame()
k = 0
for val in TE.iterrows():
    # display(val[0:2])
    # print(val[1][0])
    gene_name = val[0]
    # print([val[1][0], val[1][1], val[1][2]], [val[1][6], val[1][7], val[1][8]])
    statD4, p_valD4 = stats.ttest_rel([val[1][0], val[1][1], val[1][2]], [val[1][3], val[1][4], val[1][5]])
    statD8, p_valD8 = stats.ttest_rel([val[1][0], val[1][1], val[1][2]], [val[1][6], val[1][7], val[1][8]])
    log_pval_D4 = math.log10(p_valD4)*(-1)
    log_pval_D8 = math.log10(p_valD8)*(-1)
    
    p[gene_name] = [statD4, p_valD4, log_pval_D4, statD8, p_valD8, log_pval_D8]
p = p.T
p.columns = ['statD4', 'p_valD4', 'log_pval_D4', 'statD8', 'p_valD8', 'log_pval_D8']
display(p)


graph = pd.concat([TE_FC.iloc[:,6:9], p],axis=1)
display(graph)

D4_in_graph = graph[(graph['log_pval_D4'] > 1.301) & (graph['meanD4'] >= 0)]
D8_in_graph = graph[(graph['log_pval_D8'] > 1.301) & (graph['meanD8'] >= 0)]
D4_de_graph = graph[(graph['log_pval_D4'] > 1.301) & (graph['meanD4'] < 0)]
D8_de_graph = graph[(graph['log_pval_D8'] > 1.301) & (graph['meanD8'] < 0)]

# 원하는 gene name의 친구들만 highlight
genelist = ['Mief1', 'Ucp2', 'Crls1', 'Slc25a44', 'Pcdhgc3']
wewant = pd.DataFrame()
for val in graph.iterrows():
    for name in genelist:
        if val[0] == name:
            print(val[0], val[1][0], name)
            print(val, type(val))
            wewant[name] = val[1:9][0]

wewant = wewant.T
display(wewant)


fig = plt.figure(figsize=(6,3)) # 36, 12

gs = gridspec.GridSpec(nrows=1, ncols=2) # 1개

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])



in_color = 'red'
de_color = 'blue'
edgecolor = 'black'
wewantcolor = 'yellow'
dot_size = 40
ax1.scatter(graph['meanD4'], graph['log_pval_D4'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(D4_in_graph['meanD4'], D4_in_graph['log_pval_D4'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(D4_de_graph['meanD4'], D4_de_graph['log_pval_D4'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax1.scatter(wewant['meanD4'], wewant['log_pval_D4'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )


ax1.text(2.5, 2.35, f'{len(D4_in_graph)}', fontdict={'color':'red'})
ax1.text(-2.75, 2.35, f'{len(D4_de_graph)}', fontdict={'color':'blue'})

ax2.scatter(graph['meanD8'], graph['log_pval_D8'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_in_graph['meanD8'], D8_in_graph['log_pval_D8'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_de_graph['meanD8'], D8_de_graph['log_pval_D8'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax2.scatter(wewant['meanD8'], wewant['log_pval_D8'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax2.text(2.5, 2.95, f'{len(D8_in_graph)}', fontdict={'color':'red'})
ax2.text(-2.8, 2.95, f'{len(D8_de_graph)}', fontdict={'color':'blue'})



ax1.set_ylabel('Statistical value (-log10 P-value)')
ax1.set_xlabel('log2 TE RPKM FC')
ax1.set_title('D4 vs D0')

ax2.set_ylabel('Statistical value (-log10 P-value)')
ax2.set_xlabel('log2 TE RPKM FC')
ax2.set_title('D8 vs D0')

ax1.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax1.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')

ax2.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax2.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')


ax1.set_xlim(-3, 3)
ax2.set_xlim(-3, 3)



plt.tight_layout()
plt.savefig('Global_TE_volcano_plot.pdf', dpi=300)
plt.show()

display(wewant.iloc[:, [0, 1, 4, 7]]) # 0, 1, 4, 7

fig = plt.figure(figsize=(6,3)) # 36, 12

gs = gridspec.GridSpec(nrows=1, ncols=2) # 1개

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])



in_color = 'red'
de_color = 'blue'
edgecolor = 'black'
wewantcolor = 'red'
dot_size = 40
ax1.scatter(graph['meanD4'], graph['log_pval_D4'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax1.scatter(D4_in_graph['meanD4'], D4_in_graph['log_pval_D4'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax1.scatter(D4_de_graph['meanD4'], D4_de_graph['log_pval_D4'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(wewant['meanD4'], wewant['log_pval_D4'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax2.scatter(graph['meanD8'], graph['log_pval_D8'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax2.scatter(D8_in_graph['meanD8'], D8_in_graph['log_pval_D8'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax2.scatter(D8_de_graph['meanD8'], D8_de_graph['log_pval_D8'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(wewant['meanD8'], wewant['log_pval_D8'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax1.set_ylabel('Statistical value (-log10 P-value)')
ax1.set_xlabel('log2 TE RPKM FC')
ax1.set_title('D4 vs D0')

ax2.set_ylabel('Statistical value (-log10 P-value)')
ax2.set_xlabel('log2 TE RPKM FC')
ax2.set_title('D8 vs D0')

ax1.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax1.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')

ax2.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax2.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')


ax1.set_xlim(-3, 3)
ax2.set_xlim(-3, 3)



plt.tight_layout()
plt.savefig('Highlight_TE_volcano_plot.pdf', dpi=300)
plt.show()

display(wewant.iloc[:, [0, 1, 4, 7]]) # 0, 1, 4, 7
