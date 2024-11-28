import pandas as pd

for j in 0, 4, 8:

    com = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{j}.ribocode.riborf.csv', sep='\t')
    #uORF (uORF)
    com[lambda x: com['ORF_type'] == 'uORF'].to_csv(f'/Data_3/Suna/ORF/commonORF/type/D{j}/uORF.csv', sep='\t', index=False) 


    #dORF (dORF)
    com[lambda x: com['ORF_type'] == 'dORF'].to_csv(f'/Data_3/Suna/ORF/commonORF/type/D{j}/dORF.csv', sep='\t', index=False) 


    #ncRNA (novel : noncoding gene/transcript)
    com[lambda x: (com['ORF_type'] == 'novel')].to_csv(f'/Data_3/Suna/ORF/commonORF/type/D{j}/ncORF.csv', sep='\t', index=False) 

