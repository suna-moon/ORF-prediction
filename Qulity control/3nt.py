# %%
# 3nt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

tri = ['a', 'b', 'c']
time = ['0', '4', '8']


# %%
def threent (df):
    data = {0 : [0], 1 : [0], 2 : [0]}
    dfCount = pd.DataFrame(data)
    for i in df.loc[11:, :].index:
        if i % 3 == 0:
            dfCount.at[0, 0] += df.at[i, '0']
        elif i % 3 ==1 :
            dfCount.at[0, 1] += df.at[i, '0']
        elif i % 3 == 2 :
            dfCount.at[0, 2] += df.at[i, '0']
        

    return dfCount

k = 0
fig = plt.figure(figsize=[12, 16])

for i in time:
    for j in tri:
        k += 1

        Fsize = 18
        titleFsize = 25
        df = pd.read_csv(f'/Data_3/Suna/Adipocytes/start/D{i}{j}.csv')
        # display(df)
        # print(df['0'])
        # print(df['Unnamed: 0'])

        # display(df['0'])
        df = df.set_index(keys=['Unnamed: 0'])
        # dfStart = pd.DataFrame({0:df['0']}, index = df['Unnamed: 0']) # 열 바꾸기

        # dis/ay(df)
        dfCount = threent(df)
        # # display(dfCount)
        # # display(dfCount.T)
        # # print(dfCount.T.sum(axis=0))
        # display(dfCount)
        dfCount = dfCount.T/dfCount.T.sum(axis=0)*100
        
        ax = fig.add_subplot(3, 3, k)

        ax.bar(dfCount.index, dfCount[0], color = ['lightcoral', 'silver', 'grey'], zorder = 1000, edgecolor = 'black',  linewidth = 0.8)

        ax.set_title(f'D{i}{j}',  fontdict={'fontsize':titleFsize})
        ax.set_xticks(np.arange(0, 3, 1))
        ax.set_xticklabels(np.arange(0, 3, 1), fontdict={'fontsize':Fsize})
        ax.set_yticks(np.arange(0, 80, 10))
        ax.set_ylim(0, 70)
        ax.set_xlabel('Frame', fontdict={'fontsize':Fsize})
        ax.set_ylabel('Percentage(%)',  fontdict={'fontsize':Fsize})

        ax.grid(True, zorder = 1, linestyle='--')
        
plt.tight_layout()
plt.savefig("/Data_3/Suna/Adipocytes/fig/3nt.pdf", dpi = 300, bbox_inches = "tight")
plt.show()





