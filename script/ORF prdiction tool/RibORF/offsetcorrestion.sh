# /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/offsetCorrect.pl
# /Data_3/Suna/RibORFfile/output/D0a.sorted.sam
# /Data_3/Suna/RibORFfile/output/D0b.offsetCorrection.txt
# perl offsetCorrect.pl -r SRR1802146.mapping.sam -p offset.correction.parameters.txt -o corrected.SRR1802146.mapping.sam


# /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/readDist.pl
# /output/D${i}${j}.corrected.adipocytes.mapping.sam
# /Data_1/References/mm39_vM27/Mouse.genePred.txt
# /Data_3/Suna/RibORFfile/output/correctedReadDist/D${i}${j}
# perl readDist.pl -f corrected.SRR1802146.mapping.sam -g gencode.v28.annotation.genePred.txt -o outputDir -d 1



for i in 0 4 8
do
for j in a b c
do
perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/offsetCorrect.pl -r /Data_3/Suna/RibORFfile/output/D${i}${j}.sorted.sam -p /Data_3/Suna/RibORFfile/output/correctedOffsetData/D${i}${j}.txt -o /Data_3/Suna/RibORFfile/output/corretedSamFile/D${i}${j}.corrected.sam
perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/readDist.pl -f /Data_3/Suna/RibORFfile/output/corretedSamFile/D${i}${j}.corrected.sam -g /Data_1/References/mm39_vM27/Mouse.genePred.txt -o /Data_3/Suna/RibORFfile/output/correctedReadDist/D${i}${j} -d 1
done
done
