# %%
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

time = [0, 4, 8]
tri = ['a', 'b', 'c']

for i in time:
    for j in tri:

        a = 0
        b = 0
                
        dfStart = pd.read_csv(f'/Data_2/Jun/Adipocytes/rpf/novaseq/D{i}{j}.counts-around-start.txt', sep='\t').set_index('pos').loc[-23 : 10, '23' : '36'].T
        dfStop = pd.read_csv(f'/Data_2/Jun/Adipocytes/rpf/novaseq/D{i}{j}.counts-around-stop.txt', sep='\t').set_index('pos').loc[-33 : -3, '23' : '36'].T  # +3 으로 -30 -3 되도록
        
        a = dfStart.sum().sum() + dfStop.sum().sum()
        dfStop.columns = dfStop.columns + 3  # stop codon + 3
        dfStart = dfStart/int(a)*10000
        dfStop = dfStop/int(a)*10000


        # 그래프 그리기

        fig = plt.figure(figsize=[16, 4])

        ax1 = fig.add_subplot(1, 2, 1)
        ax1 = sns.heatmap(dfStart, cmap='hot', cbar=False, vmin=0, vmax=220)
        plt.xticks(rotation = 45)
        plt.xlabel('Position From Start Codon', fontsize=15)
        plt.ylabel('Read Length', fontsize=15)
        plt.yticks(rotation = 0)
        ax1.set_title(f'D{i}{j} Read Length Metagene (Start)', fontsize=15)

        ax2 = fig.add_subplot(1, 2, 2)
        ax2 = sns.heatmap(dfStop, cmap='hot', cbar=False, vmin=0, vmax=220)
        ax2.set_yticks([], [])
        plt.xticks(rotation = 45)
        plt.xlabel('Position From Stop Codon', fontsize=15)
        # plt.xlim([-30, 0])
        ax2.set_title(f'D{i}{j} Read Length Metagene (Stop)', fontsize=15)

        plt.savefig(f"/Data_3/Suna/Adipocytes/fig/readlengthmeatagene/D{i}{j}.pdf", dpi = 300, bbox_inches = "tight")
        plt.tight_layout()







