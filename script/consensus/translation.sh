# uORF : Ucp2 Ppp2r5a Crls1 Btg1 Slc25a44 Ppp1r15a Mief1 Vgll3 Pcdhgc3 Mgat2
# dORF : Col12a1 Dph1 Peg10 Ppp1r12a Smim10l1
# lncRNA : C4a
# conda emboss


for i in C4a_4 C4a_6 C4a_8
do
    transeq /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/finalinput/lncRNA/${i}.fasta /Data_3/Suna/ORF/Consensus/protein_seq/lncRNA/${i}.fa
done
# transeq /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/MSA_clustalOmega/uORF/Ucp2_14.fasta Ucp2.fa
# transeq /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/finalinput/uORF/Abca1.fasta
