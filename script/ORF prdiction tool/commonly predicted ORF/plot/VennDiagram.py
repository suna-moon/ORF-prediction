#ribocode : /Data_3/Suna/ORF/RiboCode/commonORF/D0.length.txt
#ribORF : /Data_3/Suna/ORF/RibORF/commonORF/D0.length.txt
#common : /Data_3/Suna/ORF/commonORF/D0.ribocode.riborf.csv

from matplotlib_venn import venn2
import pandas as pd
from matplotlib import pyplot as plt
fig, axes = plt.subplots(1, 3, figsize=[12, 4])
j = 0
for i in 0, 4, 8:
    
    codeNumber = len(pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.repre.txt'))
    orfNumber = len(pd.read_csv(f'/Data_3/Suna/ORF/RibORF/commonORF/D{i}.length.repre.txt'))
    common = len(pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{i}.ribocode.riborf.csv'))
    ax = axes[j]
    ax.set_title(f'D{i}')
    # venn2(subset = (왼쪽원-교개수, 오른쪽원-교개수, 교집합개수))
    # axes.title(f'D{i}')
    # plt.title(f'D{i}')
    j += 1
    print(codeNumber, common, orfNumber)
    venn2(subsets = (codeNumber-common, orfNumber-common, common), set_labels = ('RiboCode', 'RibORF'), ax=ax)
    
plt.savefig(f"/Data_3/Suna/ORF/commonORF/VenDiagram/VennDiagram.pdf", dpi = 300, bbox_inches = "tight")

plt.show()
