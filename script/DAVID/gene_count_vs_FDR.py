import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
raw = pd.read_csv('/Data_3/Suna/figure/4.David/output/GOTERM_ALL.txt', sep='\t', header=0)
raw['-log10FDR'] = np.log10(raw['FDR'])*(-1)
# raw
raw['Term'] = [i.split('~')[1] for i in raw['Term']]
raw.sort_values(by='Count', ascending=False)
membrane = ['organelle membrane', 'intracellular membrane-bounded organelle', 'endomembrane system', 'membrane-bounded organelle']
highlight = pd.DataFrame()
for i in raw.iterrows():
    print(i[1][1])
    for j in membrane:
        if i[1][1] == j:
            highlight[f'{i[0]}'] = raw.iloc[i[0]]
highlight = highlight.T
display(highlight)

pos_FDR = raw[raw['-log10FDR'] > 1.301]
display(pos_FDR)
# color, 'Reds'
fig = plt.figure(figsize=(4,4))

plt.scatter(raw['Count'], raw['-log10FDR'], s=raw['Count'], c=raw['-log10FDR'], cmap = 'Reds', vmin=0, vmax=2) #color = 'lightgrey') # lightgrey
# cbar = vmin = 
plt.scatter(pos_FDR['Count'], pos_FDR['-log10FDR'], s=pos_FDR['Count'], c=pos_FDR['-log10FDR'], cmap = 'Reds', linewidth=0.3, vmin=0, vmax=2) 
# size = [i for i in highlight['Count']]
# plt.scatter(highlight['Count'], highlight['-log10FDR'], s=size, color='r')



plt.axhline(y=1.301, color='black', linestyle='--', linewidth=0.7) # grey, red

plt.xlabel('Gene count')
plt.ylabel('-log10(FDR)')

plt.title('David GOTERM')
plt.tight_layout()
plt.savefig('Dvaid GOTERM.pdf', dpi=300)
plt.show()
display(pos_FDR.sort_values(by='Count', ascending=False))
