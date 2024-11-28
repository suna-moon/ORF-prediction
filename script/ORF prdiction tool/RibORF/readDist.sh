for date in 0 4 8
do
for tri in a b c
do
echo /Data_2/Jun/Adipocytes/alignments/novaseq/D${date}${tri}.sorted.bam
# sam file로 변환
# samtools view /Data_2/Jun/Adipocytes/alignments/novaseq/D${date}${tri}.sorted.bam > /Data_3/Suna/RibORFfile/output/D${date}${tri}.sorted.sam

# readDist.pl 실행
perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/readDist.pl -f /Data_3/Suna/RibORFfile/output/D${date}${tri}.sorted.sam -g /Data_1/References/mm39_vM27/Mouse.genePred.txt -o output/D${date}${tri} -d 28,29,30,31,32,33,34 -l 30 -r 50
echo ${date}${tri} finish
# default
done
done
# perl /Data_3/Suna/RibORFfile/RibORF/RibORF.2.0/readDist.pl -f D0a.sorted.sam -g /Data_1/References/mm39_vM27/Mouse.genePred.txt -o output -d 28,29,30,31,32,33,34 -1 30 -r 50

