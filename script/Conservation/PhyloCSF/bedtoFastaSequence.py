import pandas as pd
cate = ['uORF', 'dORF', 'lncRNA']
for i in cate: # human cordinate만 가져오기 col7-12
    df = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/compare/dup/{i}_human_cordinate_duplicate.txt', sep='\t', header=None)
    df.iloc[:, 6:12].to_csv(f'/Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/bed6_Liftover/human/{i}.bed', sep='\t', header=None, index=None)
    df.iloc[:, 0:6].to_csv(f'/Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/bed6_Liftover/mouse/{i}.bed', sep='\t', header=None, index=None)
    display(df)
