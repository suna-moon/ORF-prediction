from pymsaviz import MsaViz
from Bio.Seq import Seq
import inspect
print(inspect.getfile(MsaViz))
top = ['Ucp2', 'Ppp2r5a', 'Crls1', 'Btg1', 'Slc25a44', 'Ppp1r15a', 'Mief1']
bottom = ['Vgll3', 'Pcdhgc3', 'Mgat2']

dORF_gene_list = ['Col12a1', 'Dph1', 'Peg10', 'Ppp1r12a', 'Smim10l1']
lncRNA_list = ['C4a_4', 'C4a_6', 'C4a_8']

# /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/finalinput/uORF
gene_list = lncRNA_list
loc = 'lncRNA'
for i in gene_list:
    msa_file = f'/Data_3/Suna/ORF/Consensus/MSA_align/{loc}/{i}.fa'
    mv = MsaViz(msa_file, wrap_length=60, show_consensus=True, color_scheme='Clustal') # , show_grid=True
    mv.savefig(f'/Data_3/Suna/ORF/Consensus/output/{loc}/{i}.pdf')
