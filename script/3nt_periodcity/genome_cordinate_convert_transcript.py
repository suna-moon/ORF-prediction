import pandas as pd
info = pd.read_csv('/Data_3/Suna/ORF/featureCounts/uORF(ORF_pos_norep).txt', sep='\t', header=0)
orf_cordinate = pd.read_csv('/Data_3/Suna/ORF/featureCounts/uORFcordinate.txt', sep='\t', header=0)

display(orf_cordinate)
# transcript_id 정보
# merge 용 ORFstart
newORFstartForMerge = [int(i.split('.')[0]) for i in orf_cordinate["ORFstart"]]
orf_cordinate['ORF_gstart'] = newORFstartForMerge
info = info.iloc[:, [3, 5, 6, 10, 13, 14]]

merge = pd.merge(left=info, right=orf_cordinate, how='right',
                 left_on='ORF_gstart', right_on='ORF_gstart')
display(merge)

# parsed
parsed = pd.read_csv('/Data_3/Suna/parsing/data.txt', sep='\t')
exon = pd.merge(left=merge, right=parsed.iloc[:, [0, 11, 12]], how='left',
                left_on='transcript_id', right_on='transcript_id')
for i, name in enumerate(exon.columns):
    print(i, name)
exon = exon.iloc[:, [0, 1, 2, 3, 4, 5, 8, 11, 12, 13]]

exon.to_csv('before_transcript_convert.txt', sep='\t', index=False)
display(exon)

# 계산
# print(exon.columns)
# 1~10
# ['transcript_id', 'gene_id', 'gene_name_x', 'ORF_length', 'ORF_gstart', 'ORF_gstop', 'chrom_x', 'strand_x', 'exon_start', 'exon_stop']
# ['gene_id', 'GeneID', 'transcript_id', 'ORFstart', 'ORFstop', 'strand_x', 'exon_start', 'exon_stop']
newdf = pd.DataFrame(columns=['transcript_id', 'gene_id', 'gene_name', 'ORF_length', 'ORF_gstart', 'ORF_gstop', 'chrom', 'strand', 'exon_start', 'exon_stop','t_start', 't_stop'])
for num, row in enumerate(exon.itertuples(name=None)):
    orf_start = int(row[5])
    orf_stop = int(row[6])
    exon_total = [[int(i), int(j)] for i, j in zip(row[9].split('.'), row[10].split('.'))]
    # '+' strand
    if row[8] == '+':
        start_i = [[orf_start-i, k, i] for k, (i, j) in enumerate(exon_total) if orf_start >= i if orf_start <= j] # 간격, n-1 번째 exon에 start 위치
        stop_i = [[orf_stop-i, k, i] for k, (i, j) in enumerate(exon_total) if orf_stop >= i if orf_stop <= j] # 간격, n-1 번째 exon에 start 위치

        # 계산
        # transcript cordinate
        t_cordinate = []
        for k, (i, j) in enumerate(exon_total):
            if k == 0:
                t_cordinate.append([1, j-i+1])
                next = j-i+1
            elif k == stop_i[0][1] + 1:  # 전체 exon의 transcript cordinate가 필요하면 elif 부분 삭제
                break
            else:
                t_cordinate.append([next+1, next+1+j-i]) # (next+1) + 간격
                next = next+1+j-i

        # orf cordinate 
        start = t_cordinate[start_i[0][1]][0] + start_i[0][0]
        stop = t_cordinate[stop_i[0][1]][0] + stop_i[0][0]
        row_list = list(row[1:11])
        row_list.extend([start, stop])
        newdf.loc[len(newdf)] = row_list

    

    # '-' strand
    if row[8] == '-':
        print(row)
        orf_start = int(row[5])
        orf_stop = int(row[6])
        start_i = [[j-orf_start, k, i, orf_start] for k, (i, j) in enumerate(exon_total) if orf_start >= i if orf_start <= j] # 간격, n-1 번째 exon에 start 위치
        stop_i = [[j-orf_stop, k, i, orf_stop] for k, (i, j) in enumerate(exon_total) if orf_stop >= i if orf_stop <= j] # 간격, n-1 번째 exon에 start 위치


        # 계산
        # transcript cordinate
        t_cordinate = []
        for k, (i, j) in enumerate(exon_total):
            if k == 0:
                print(i, j, i-j+1)
                t_cordinate.append([1, j-i+1])
                next = j-i+1
            elif k == stop_i[0][1] + 1:  # 전체 exon의 transcript cordinate가 필요하면 elif 부분 삭제
                break
            else:
                t_cordinate.append([next+1, next+1+j-i]) # (next+1) + 간격
                next = next+1+j-i
        # orf cordinate 
        start = t_cordinate[start_i[0][1]][0] + start_i[0][0]
        stop = t_cordinate[stop_i[0][1]][0] + stop_i[0][0]
        row_list = list(row[1:11])
        row_list.extend([start, stop])
        newdf.loc[len(newdf)] = row_list
display(newdf)
newdf.to_csv('transcript_convert.txt', sep='\t', index=False)
