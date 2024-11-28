# dORF : Col12a1 Dph1 Peg10 Ppp1r12a Smim10l1
# uORF : Ucp2 Ppp2r5a Crls1 Btg1 Slc25a44 Ppp1r15a Mief1 Vgll3 Pcdhgc3 Mgat2
# lncRNA : C4a_test
# codna gencode

for filename in C4a_4 C4a_6 C4a_8
    do
        clustalo -i /Data_3/Suna/ORF/Consensus/protein_seq/lncRNA/${filename}.fa -o /Data_3/Suna/ORF/Consensus/MSA_align/lncRNA/${filename}.fa --outfmt=fa --force
    done
