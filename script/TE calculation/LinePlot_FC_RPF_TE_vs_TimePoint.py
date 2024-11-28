# Data
# uORF RPKM log2 FC /Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_uORF_log2_FC_RPF.txt
# CDS RPKM log2 /Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_CDS_log2_FC_RPF.txt
# RNA RPKM 1og2 /Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_log2_FC_RNA.txt

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import statistics as st

uORF = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_uORF_log2_FC_RPF.txt', sep='\t', header=0, index_col='Geneid')
CDS = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_CDS_log2_FC_RPF.txt', sep='\t', header=0, index_col='Geneid')
RNA = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/TE/RPKM_log2_FC_RNA.txt', sep='\t', header=0, index_col='gene_name')

# display(uORF, CDS, RNA)
uORF_want = pd.DataFrame()
CDS_want = pd.DataFrame()
RNA_want = pd.DataFrame()
genelist = ['Mief1', 'Ucp2', 'Slc25a44', 'Pcdhgc3', 'Crls1']
# print(CDS.loc['Mief1'])
for i, name in enumerate(genelist):
    uORF_want[name] = uORF.loc[f'{name}']
    CDS_want[name] = CDS.loc[f'{name}']
    RNA_want[name] = RNA.loc[f'{name}']

# display(uORF)
uORF = [uORF.loc[name] for name in genelist]

uORF_want = uORF_want.T.iloc[:, 18:24]
# uORF_want = uORF_want.T.iloc[:, list(range(9, 12)) + list(range(18, 24))]
CDS_want = CDS_want.T.iloc[:, 18:24]
RNA_want = RNA_want.T.iloc[:, 18:24]

display(uORF_want, CDS_want, RNA_want)

# TE caculation
TE_CDS = CDS_want-RNA_want
TE_uORF =uORF_want-RNA_want
TE_uORF['meanD0'] = [0, 0, 0, 0, 0]
TE_CDS['meanD0'] = [0, 0, 0, 0, 0]
uORF_want['meanD0'] = [0, 0, 0, 0, 0]
CDS_want['meanD0'] = [0, 0, 0, 0, 0]

# mean
for i, j in enumerate([4, 8]):
    TE_CDS[f'meanD{j}'] = TE_CDS.iloc[:, i*3:(i+1)*3].mean(axis=1)
    TE_uORF[f'meanD{j}'] = TE_uORF.iloc[:, i*3:(i+1)*3].mean(axis=1)
    CDS_want[f'meanD{j}'] = CDS_want.iloc[:, i*3:(i+1)*3].mean(axis=1)
    uORF_want[f'meanD{j}'] = uORF_want.iloc[:, i*3:(i+1)*3].mean(axis=1)


dfs = [TE_uORF, TE_CDS, uORF_want, CDS_want]
for df in dfs:
    for i, j in enumerate([4, 8]):
        std = []
        for val in df.iterrows():
            sample = []
            if j == 4:
                for k in range(0,3):
                    sample.append(val[1][k])
                std.append(st.stdev(sample))
            if j == 8:
                for k in range(3,6):
                    sample.append(val[1][k])
                std.append(st.stdev(sample))
        df[f'stdD{j}'] = std
        df['stdD0'] = [0, 0, 0, 0, 0]
TE_uORF, TE_CDS, uORF_want, CDS_want = dfs
    #         val[1][k]
    # print([val[1][k] for val in TE_CDS.iterrows() for k in range(0,3)])
    # TE_CDS[f'meanD{j}'] = st.stdev([])
# print(st.stdev([-2.313632904670849, -1.9030522164575512, -2.110358812526689]))

display(TE_uORF, TE_CDS, uORF_want, CDS_want)

# line chart for log2 FC (D0~8) graph

# gene_name = 'Mief1'
for gene_name in genelist:
    fig = plt.figure(figsize=(6,3))#, dpi=300) # 36, 12

    gs = gridspec.GridSpec(nrows=1, ncols=2) # 1개

    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    dot_size = 8
    line_width = 2

    barWideth = 0.9
    x = ['D0', 'D4', 'D8']
    ax1.plot(x, [i for i in uORF_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'red', linewidth=line_width)
    ax1.scatter(x, [i for i in uORF_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'red', s=dot_size)
    ax1.errorbar(x, [i for i in uORF_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'red',
                yerr=[i for i in uORF_want.loc[:, ['stdD0', 'stdD4', 'stdD8']].loc[gene_name]], capsize=3)
    ax1.plot(x, [j for j in CDS_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]], color = 'blue', linewidth=line_width)
    ax1.scatter(x, [j for j in CDS_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'blue', s=dot_size)
    ax1.errorbar(x, [i for i in CDS_want.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'blue', 
                yerr=[i for i in CDS_want.loc[:, ['stdD0', 'stdD4', 'stdD8']].loc[gene_name]], capsize=3)


    # TE
    ax2.plot(x, [i for i in TE_uORF.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]], color = 'red', linewidth=line_width)
    ax2.plot(x, [j for j in TE_CDS.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]], color = 'blue', linewidth=line_width)
    ax2.scatter(x, [i for i in TE_uORF.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]], color = 'red', s=dot_size)
    ax2.scatter(x, [j for j in TE_CDS.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]], color = 'blue', s=dot_size)
    ax2.errorbar(x, [i for i in TE_uORF.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'red', 
                yerr=[i for i in TE_uORF.loc[:, ['stdD0', 'stdD4', 'stdD8']].loc[gene_name]], capsize=3)
    ax2.errorbar(x, [i for i in TE_CDS.loc[:, ['meanD0', 'meanD4', 'meanD8']].loc[gene_name]],color = 'blue', 
                yerr=[i for i in TE_CDS.loc[:, ['stdD0', 'stdD4', 'stdD8']].loc[gene_name]], capsize=3)

    ax1.set_ylabel('log2 RPF RPKM FC')
    ax1.set_xlabel('Time')
    ax1.set_title('RPF')
    # ax1.set_xticks([1.5, 4.5])
    ax1.set_xticklabels(['D0', 'D4', 'D8'])

    ax2.set_ylabel('log2 TE RPKM FC')
    ax2.set_xlabel('Time')
    ax2.set_title('TE')
    # ax2.set_xticks([1.5, 4.5])
    ax2.set_xticklabels(['D0', 'D4', 'D8'])


    st = fig.suptitle(f'{gene_name}')

    plt.gca().legend(('uORF', 'CDS'), loc='center', bbox_to_anchor=(1.4, 0.5), edgecolor='white')
    ax1.axhline(0, color='black', linewidth=0.7, linestyle='--')
    ax2.axhline(0, color='black', linewidth=0.7, linestyle='--')

    a, b = ax1.get_ylim()
    c ,d = ax2.get_ylim()

    y_max = max(b, d)
    y_min = min(a, c)

    ax1.set_ylim(y_min, y_max)
    ax2.set_ylim(y_min, y_max)
    plt.tight_layout()
    plt.savefig(f'/Data_3/Suna/figure/5.RPF_TE_for_each_gene/plot_output/{gene_name}.pdf', dpi=300)
    plt.show()
