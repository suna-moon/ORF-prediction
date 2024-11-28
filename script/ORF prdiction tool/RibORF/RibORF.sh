# /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/ribORF.pl
# /Data_3/Suna/RibORFfile/output/D${i}${j}.corrected.adipocytes.mapping.sam


# perl ribORF.pl -f corrected.SRR1802146.mapping.sam -c candidateORF.genepred.txt -o outputDir
# perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/ribORF.pl -f /Data_3/Suna/RibORFfile/output/corretedSamFile/D0a.corrected.sam -c /Data_3/Suna/RibORFfile/output/candidateORF.genepred.txt -o result/D${i}${j}
# /Data_3/Suna/RibORFfile/output/corretedSamFile/D0a.corrected.sam

# # test file D0a 용
# perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/ribORF.pl -f /Data_3/Suna/RibORFfile/output/corretedSamFile/D0a.corrected.sam -c /Data_3/Suna/RibORFfile/output/candidateORF.genepred.txt -o result/D0a
# # for문 용
# perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/ribORF.pl -f /Data_3/Suna/RibORFfile/output/corretedSamFile/D${i}${j}.corrected.sam -c /Data_3/Suna/RibORFfile/output/candidateORF.genepred.txt -o result/D${i}${j}
for i in 0 4 8
do
for j in a b c
do
perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/ribORF.pl -f /Data_3/Suna/RibORFfile/output/corretedSamFile/D${i}${j}.corrected.sam -c /Data_3/Suna/RibORFfile/output/candidateORF.genepred.txt -o result/D${i}${j}
done
done
