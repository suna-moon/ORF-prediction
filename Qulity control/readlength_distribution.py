# %%
import matplotlib.pyplot as plt
import pysam
import pandas as pd


def count_elements(lst):
    counts = {}
    for item in lst:
        if item in counts:
            counts[item] += 1
        else:
            counts[item] = 1
    return counts

readLength = []
countList = []

tri = ['a', 'b', 'c']
time = ['0', '4', '8']

for i in time:
    for j in tri:
        filepath = f'/Data_2/Jun/Adipocytes/tr-aln/novaseq/D{i}{j}.rep.bam'
        bamfile = pysam.AlignmentFile(filepath, 'rb')
        for line in bamfile:
            readLength.append(line.reference_length)
        countList.append(count_elements(readLength))

# display(countList)

# %%
import numpy as np

fig = plt.figure(figsize=[16,12])
dict = {0:1, 1:2, 2:3, 3:4}

k = 0
Fsize = 20
xtickslabelSize = 20
titleLabel = 30

tri = ['a', 'b', 'c']
time = ['0', '4', '8']

for i in time:
    for j in tri:
        k += 1
        ax = fig.add_subplot(3, 3, k)

        x = countList[k-1].keys()
        values = list(countList[k-1].values())
        total = sum(countList[k-1].values())
        y = [val/total*10 for val in values]
        # y = [val/s*100 for val in values]

        ax.bar(x, y, zorder = 100, width=1, color = 'green', edgecolor = 'black', linewidth = 1.5 )
        ax.set_xlim(20, 40)
        ax.set_ylim(0, 3)
        ax.set_title(f'D{i}{j}', fontdict={'fontsize':titleLabel})
        ax.grid(True, zorder = 1, linestyle='--')
        ax.set_xticks(np.arange(20, 44, 5))
        ax.set_xticklabels(np.arange(20, 44, 5), fontdict={'fontsize':xtickslabelSize})
        ax.set_yticks(np.arange(0, 4, 1))
        ax.set_xlabel('Read Length',  fontdict={'fontsize':Fsize})
        ax.set_ylabel('Normalized Counts', fontdict={'fontsize':Fsize})
        
plt.tight_layout()
plt.savefig("/Data_3/Suna/Adipocytes/fig/readlength.pdf", dpi = 300, bbox_inches = "tight")
plt.show()


