import pandas as pd
for i in 0, 4, 8:
    orf = pd.read_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.txt', sep='\t')
    ref = pd.read_csv(f'/Data_1/References/mm39_vM27/representative-isoforms.txt', sep='\t', header=None).iloc[:, 1]
    merge = pd.merge(orf,ref, how='inner', left_on='transcript_id', right_on=1).drop(1, axis=1).to_csv(f'/Data_3/Suna/ORF/RiboCode/commonORF/D{i}.length.repre.txt', sep='\t', index=False)
    display(merge)

# display(orf.sort_values(by='TranscriptID',ascending=True).to_csv('/Data_3/Suna/ORF/RibORF/commonORF/isoformFilterTest.txt', sep='\t'))
