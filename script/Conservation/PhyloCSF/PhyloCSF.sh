# PHYLOCSF="/home/suna/miniconda3/envs/PhyloCSF/bin" # /home/suna/miniconda3/envs/PhyloCSF/bin
# option 없는 dafault 값

filename=`find /Data_3/Suna/ORF/PhyloCSF/InputSequenc_Fasta/MSA_clustalOmega/lncRNA -type f`
for i in ${filename}
do
PhyloCSF 30mammals_mm39 ${i} --removeRefGaps >> lncRNA.txt 
done


# PhyloCSF 30mammals_mm39 Ttc23_234.fasta >> uORF.txt 

# --strategy=fixed --frames=1 --bls --ancComp --removeRefGaps
#30mammals_mm39
# /Data_3/Suna/ORF/PhyloCSF/C4a_4.fasta
# $PHYLOCSF --strategy=fixed --frames=1 --bls --ancComp --removeRefGaps 100vertebrates /tmp/phylo.fasta
