import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec

import math
from scipy import stats

# data RPKM

# uORF RPKM data (not fold change)
rpkm = pd.read_csv('/Data_3/Suna/ORF/featureCounts/RPKM.txt', sep='\t', header =0)
# uORF RPKM FC data
fc = pd.read_csv('/Data_3/Suna/ORF/featureCounts/RPKM_FC.txt', sep='\t', header=0, index_col='Geneid')
fc = fc.iloc[:,18:25]

fc['meanD4'] = fc.iloc[:, 0:3].mean(axis=1)
fc['meanD8'] = fc.iloc[:, 3:6].mean(axis=1)
display(fc)

# paired t-test
p = pd.DataFrame()
for val in rpkm.iterrows():
    gene_name = val[1][0]
    statD4, p_valD4 = stats.ttest_rel([val[1][1], val[1][2], val[1][3]], [val[1][4], val[1][5], val[1][6]])
    statD8, p_valD8 = stats.ttest_rel([val[1][1], val[1][2], val[1][3]], [val[1][7], val[1][8], val[1][8]])
    log_pval_D4 = math.log10(p_valD4)*(-1)
    log_pval_D8 = math.log10(p_valD8)*(-1)
    p[gene_name] = [statD4, p_valD4, log_pval_D4, statD8, p_valD8, log_pval_D8]
p = p.T
p.columns = ['statD4', 'p_valD4', 'log_pval_D4', 'statD8', 'p_valD8', 'log_pval_D8']
display(p)

graph = pd.concat([fc, p],axis=1)
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


ax1.text(5, 2.6, f'{len(D4_in_graph)}', fontdict={'color':'red'})
ax1.text(-6, 2.6, f'{len(D4_de_graph)}', fontdict={'color':'blue'})

ax2.scatter(graph['meanD8'], graph['log_pval_D8'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_in_graph['meanD8'], D8_in_graph['log_pval_D8'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_de_graph['meanD8'], D8_de_graph['log_pval_D8'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
# ax2.scatter(wewant['meanD8'], wewant['log_pval_D8'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax2.text(6, 3.5, f'{len(D8_in_graph)}', fontdict={'color':'red'})
ax2.text(-7, 3.5, f'{len(D8_de_graph)}', fontdict={'color':'blue'})



ax1.set_ylabel('Statistical value (-log10 P-value)')
ax1.set_xlabel('log2 RPF RPKM FC')
ax1.set_title('D4 vs D0')

ax2.set_ylabel('Statistical value (-log10 P-value)')
ax2.set_xlabel('log2 RPF RPKM FC')
ax2.set_title('D8 vs D0')

ax1.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax1.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')

ax2.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax2.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')


ax1.set_xlim(-6.5, 6.5)
ax2.set_xlim(-7.5, 7.5)



plt.tight_layout()
plt.savefig('Global_volcano_plot.pdf', dpi=300)
plt.show()

fig = plt.figure(figsize=(6,3)) # 36, 12

gs = gridspec.GridSpec(nrows=1, ncols=2) # 1개

ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])



in_color = 'lightgrey'
de_color = 'lightgrey'
edgecolor = 'black'
wewantcolor = 'r'
dot_size = 40
ax1.scatter(graph['meanD4'], graph['log_pval_D4'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(D4_in_graph['meanD4'], D4_in_graph['log_pval_D4'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(D4_de_graph['meanD4'], D4_de_graph['log_pval_D4'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax1.scatter(wewant['meanD4'], wewant['log_pval_D4'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax2.scatter(graph['meanD8'], graph['log_pval_D8'], color='lightgrey', s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_in_graph['meanD8'], D8_in_graph['log_pval_D8'], color=in_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(D8_de_graph['meanD8'], D8_de_graph['log_pval_D8'], color=de_color, s=dot_size, edgecolor = 'black', linewidths=0.3)
ax2.scatter(wewant['meanD8'], wewant['log_pval_D8'], color=wewantcolor, s=dot_size, edgecolor = 'black', linewidths=0.3 )

ax1.set_ylabel('Statistical value (-log10 P-value)')
ax1.set_xlabel('log2 RPF RPKM FC')
ax1.set_title('D4 vs D0')

ax2.set_ylabel('Statistical value (-log10 P-value)')
ax2.set_xlabel('log2 RPF RPKM FC')
ax2.set_title('D8 vs D0')

ax1.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax1.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')

ax2.axhline(1.301, 0, 1, color = 'black', linewidth=0.7, linestyle='--')
ax2.axvline(0, 0, 1, color = 'black', linewidth=0.7, linestyle='--')


ax1.set_xlim(-6.5, 6.5)
ax2.set_xlim(-7.5, 7.5)



plt.tight_layout()
plt.savefig('Higlight_volcano_plot.pdf')
plt.show()
