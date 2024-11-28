# featureCounts [options] -a <annotation_file> -o <output_file> input_file1 [input_file2] ...
# gene ?, transcript ????

# featureCounts -a /ORF/featureCounts/uORFcordinate.saf -o featureCountes_output.txt bam1.file bam1.file ...
# # option
# -T 4 #  cpu 사용 지정
# -s 1 # 1 for stranded, 2 for reversely stranded, and 0 for unstranded., default = 0
# -t string # annotation이 GTF file일 경우 GTF의 feature 중 선택 exon, transcript, UTR, CDS , transcript/exon default = exon
# -p # paired end의 경우 필요 X
# -M # multiple location에 mapping되는 read(sequence similarities)를 포함 optional  X
# -O # multiple overlapping (alternative splicing, transcript isoform) optional x
# -g # GTF에서만 사용
# -F SAF gtf/saf fileformat defacult = GTF
# -o # output file format 
# default = sam, bam 사용은 옵션 필요X, transcript align파일 사용 => rep.bam

# filename=`find /Data_2/Jun/Adipocytes/alignments/novaseq -name "*.sorted.bam"`+" "
# echo $filename

# /Data_2/Jun/Adipocytes/alignments/novaseq/D0a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D0b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D0c.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4c.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8c.sorted.bam
# /Data_2/Jun/Adipocytes/alignments/novaseq/D0a.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D0b.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D0c.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D4a.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D4b.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D4c.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D8a.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D8b.sorted.bam 
# /Data_2/Jun/Adipocytes/alignments/novaseq/D8c.sorted.bam

featureCounts -T 4 -s 1 -F SAF -a /Data_3/Suna/ORF/featureCounts/lncRNAcordinate.saf -o featureCountes_output_lncRNA.txt /Data_2/Jun/Adipocytes/alignments/novaseq/D0a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D0b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D0c.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D4c.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8a.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8b.sorted.bam /Data_2/Jun/Adipocytes/alignments/novaseq/D8c.sorted.bam
