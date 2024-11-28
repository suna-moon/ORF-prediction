#fasta : 
    # mouse : /Data_1/References/mm39_vM27/Mouse.fasta
    # human : /Data_1/References/hg38_v42/Human.fasta
#bed : /Data_3/Suna/ORF/PhyloCSF/InputSequence(Fasta)/bed6(Liftover)/mouse/uORF.bed


for j in uORF dORF lncRNA
do

    bedtools getfasta -fi /Data_1/References/mm39_vM27/Mouse.fasta -bed /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/bed6_Liftover/mouse/${j}.bed -fo /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/mouse/${j}.fasta -name -s
    bedtools getfasta -fi /Data_1/References/hg38_v42/Human.fasta -bed /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/bed6_Liftover/human/${j}.bed -fo /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/human/${j}.fasta -name -s

done
