# %%
# metagene plot

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tri = ['a', 'b', 'c']
time = ['0', '4', '8']

for i in time:
    for j in tri:

        Fsize = 15

        dfStart = pd.read_csv(f'/Data_3/Suna/Adipocytes/start/D{i}{j}.csv')
        dfStop = pd.read_csv(f'/Data_3/Suna/Adipocytes/stop/D{i}{j}.csv')
        # display(dfStart)
        # print(dfStart.index)
        # print(dfStart['Unnamed: 0'])
        # print(dfStart['0'])
        
        a = dfStart.sum(axis=0) + dfStop.sum(axis=0)
        dfStart['0'] = dfStart['0']/a[1]*100  # Normalization
        dfStop['0'] = dfStop['0']/a[1]*100

        fig = plt.figure(figsize=(8,3))

        ax1 = fig.add_subplot(1, 2, 1)
        plt.grid(True)  
        plt.xlabel('Relative Position From Start Codon', fontdict={'fontsize':Fsize})
        plt.ylabel('Normalized Counts', fontdict={'fontsize':Fsize})

        ax2 = fig.add_subplot(1, 2, 2)
        plt.grid(True)  
        plt.xlabel('Relative Position From Stop Codon', fontdict={'fontsize':Fsize})
        plt.suptitle(f'D{i}{j}', fontsize=25)  # 전체 title # font size
    

        ax1.plot(dfStart['Unnamed: 0'], dfStart['0'])
        ax2.plot(dfStop['Unnamed: 0'], dfStop['0'])
        ax1.set(xlim = [-100, 100])
        ax2.set(xlim = [-100, 100])
        ax1.set(ylim = [0, 8])
        ax2.set(ylim = [0, 8])
        plt.tight_layout()  # 간격 조절

        plt.savefig(f"/Data_3/Suna/Adipocytes/fig/D{i}{j}.pdf", dpi = 300, bbox_inches = "tight")
        plt.show()



