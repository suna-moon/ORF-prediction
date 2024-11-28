# metaplots -a RiboCode_annot -r /Data_2/Jun/Adipocytes/tr-aln/novaseq/D0a.rep.bam

for i in 0 4 8
do
for j in a b c
do
metaplots -a RiboCode_annot -r /Data_2/Jun/Adipocytes/tr-aln/novaseq/D${i}${j}.rep.bam -o /Data_3/Suna/RiboCode/output/metaplot/D${i}${j}
done
done
