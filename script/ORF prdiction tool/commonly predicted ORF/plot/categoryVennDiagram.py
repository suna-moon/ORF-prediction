#ribocode : /Data_3/Suna/ORF/RiboCode/commonORF/D0.length.txt #ORF_type : uORF, dORF, novel
#ribORF : /Data_3/Suna/ORF/RibORF/commonORF/D0.length.txt # candidateORFType : uORF, dORF, noncoding
#common : /Data_3/Suna/ORF/commonORF/D0.ribocode.riborf.csv # /Data_3/Suna/ORF/commonORF/type/D0/dORF.csv 여기도 있음 filter 할 필요 없이 그냥 len만 주면 됨!
from matplotlib_venn import venn2
import pandas as pd
from matplotlib import pyplot as plt
# transcript_type의 lncRNA
day = [0, 4, 8]

for code, orf in zip(['uORF', 'dORF', 'novel'], ['uORF', 'dORF', 'noncoding']):
    fig = plt.figure(figsize=[10, 3])
    k = 0

    for j, i in enumerate(day):
        ### RiboCode
        codeData = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.repre.txt', sep='\t').loc[:, ['ORF_type']]
        if orf == 'noncoding':
            codeData = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.repre.txt', sep='\t').loc[:, ['transcript_type']]
            codeNumber = len(codeData[codeData['transcript_type'] == 'lncRNA'])
        else:
            codeData = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.repre.txt', sep='\t').loc[:, ['ORF_type']]
            codeNumber = len(codeData[codeData['ORF_type'] == f'{code}'])
        
        ### RibORF  # RibORf file에 transcript type 정보 X => parsing file에서 가져옴(using pd.merge())
        if orf == 'noncoding': 
            orfData = pd.read_csv(f'/Data_3/Suna/ORF/RibORF/commonORF/D{i}.length.repre.txt', sep='\t').loc[:, ['TranscriptID', 'candidateORFType']]
            orfData = orfData[orfData['candidateORFType'] == f'{orf}']
            parsing = pd.read_csv(f'/Data_3/Suna/parsing/data.txt', sep='\t').loc[:, ['transcript_id', 'gene_type']]
            com1 = pd.merge(left=orfData, right=parsing,
                   how='inner',
                   left_on=['TranscriptID'],
                   right_on=['transcript_id']
                   )
            commonNumber = com1[com1['gene_type'] == 'lncRNA']
        else:
            orfData = pd.read_csv(f'/Data_3/Suna/ORF/RibORF/commonORF/D{i}.length.repre.txt', sep='\t').loc[:, ['candidateORFType']]
            orfNumber = len(orfData[orfData['candidateORFType'] == f'{orf}'])

        ### Common
        if orf == 'noncoding':
            commonData = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/type/D{i}/ncORF.csv', sep='\t').loc[:, ['transcript_type']]
            commonNumber = len(commonData[commonData['transcript_type'] == 'lncRNA'])
        else:
            commonData = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{i}.ribocode.riborf.csv', sep='\t').loc[:, ['ORF_type']]
            commonNumber = len(commonData[commonData['ORF_type'] == f'{code}'])


        ax = fig.add_subplot(1, 3, k+1)
        k += 1
        titleFsize=15
        ax.set_title(f'D{i}', fontdict={'fontsize':titleFsize})
        venn2(subsets = (codeNumber-commonNumber, orfNumber-commonNumber, commonNumber), set_labels = ('RiboCode', 'RibORF'), ax=ax)
    suptitleFsize = 20
    if orf == 'noncoding':
        plt.suptitle('lncRNA', fontsize=suptitleFsize)
    else:
        plt.suptitle(f'{orf}', fontsize=suptitleFsize)  
    plt.savefig(f"/Data_3/Suna/ORF/commonORF/categoryVennDiagram{code}.pdf", dpi = 300, bbox_inches = "tight")
    plt.show()
