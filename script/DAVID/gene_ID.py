import pandas as pd

# /Data_3/Suna/ORF/featureCounts/CDS_RPF/name_id.txt
# /Data_3/Suna/ORF/featureCounts/CDS_RPF/name_id_dORF.txt

uorf = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/name_id.txt', sep='\t', header=0)
dorf = pd.read_csv('/Data_3/Suna/ORF/featureCounts/CDS_RPF/name_id_dORF.txt', sep='\t', header=0)

uorf['del_version'] = [i[1][0].split('.')[0] for i in uorf.iterrows()]
dorf['del_version'] = [i[1][0].split('.')[0] for i in dorf.iterrows()]

uorf['del_version'].to_csv('uORF_gene_id.txt', header=None, index=None)
dorf['del_version'].to_csv('dORF_gene_id', header=None, index=None)
display(uorf, dorf)
