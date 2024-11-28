# 파씽파일에서 exon 위치 정보 불러오기

# mergeORF(day)file의 Col 4: transcript_id
# 위치 : /Data_3/Suna/ORF/mergeORF(day)/dORF.csv
# parsing file의 exon_start, exon_end 정보 Col 12: exon_start, Col 13: exon_stop
# 위치 : /Data_3/Suna/parsing/data.txt

import pandas as pd

cate = ['uORF', 'dORF', 'lncRNA']

parsing = pd.read_csv('/Data_3/Suna/parsing/data.txt', sep='\t').iloc[:, [0, 11, 12]]
for i in cate:
    mergeORF = pd.read_csv(f'/Data_3/Suna/ORF/mergeORF(day)/{i}.csv', sep='\t')
    writeORF = pd.merge(left=mergeORF, right=parsing,
                        how = 'left',
                        left_on=['transcript_id'],
                        right_on=['transcript_id']
    )
    writeORF.to_csv(f'/Data_3/Suna/ORF/mergeORF(day)/{i}(exon).csv', sep='\t', index=False)
    display(writeORF)
  

######## ORF seq만 골라내기

# 위치 /Data_3/Suna/ORF/mergeORF(day)/dORF(exon).csv
# Col 14: ORF_gstart Col 15: ORF_gstop, Col 52: exon_start, Col 53: exon_stop
# Col 10: strand_x (+/-)

cate = ['uORF', 'dORF', 'lncRNA']

for i in cate:
    # position = pd.read_csv(f'/Data_3/Suna/ORF/mergeORF(day)/{i}(exon).csv', sep='\t')
    with open(f'/Data_3/Suna/ORF/mergeORF(day)/{i}(exon).csv', 'r') as f, open(f'/Data_3/Suna/ORF/mergeORF(day)/{i}(ORF_pos).csv', 'w') as out:

        for line in f:
            # print(line)
            ORFstartPosition = []
            ORFstopPosition = []
            rawORF = line.rstrip('\n').split('\t')

            if rawORF[0] == 'ORF_ID':
                rawORF.extend(['ORFstart', 'ORFstop'])
                print(rawORF)
                out.write('\t'.join(map(str, rawORF)) + '\n')


            if rawORF[9] == '+':
                exonstart = rawORF[51].split('.')
                exonstop = rawORF[52].split('.')
                startORF = int(rawORF[13])
                stopORF = int(rawORF[14])

                position = [int(exon) for exon in exonstart+exonstop if int(exon) <= stopORF and int(exon) >= startORF]
                position.extend([startORF, stopORF])
                position.sort()  # startORF, exonstop, exonstart, stopORF

                for num, pos in enumerate(position):
                    if num%2 == 0:
                        ORFstartPosition.append(pos) # 짝수 : start
                    elif num%2 == 1:
                        ORFstopPosition.append(pos) # 홀수 : stop
                ORFstartPosition = '.'.join(map(str, ORFstartPosition))
                ORFstopPosition ='.'.join(map(str, ORFstopPosition))

                rawORF.extend([ORFstartPosition, ORFstopPosition])
                out.write('\t'.join(map(str, rawORF)) + '\n')

            elif rawORF[9] == '-':
                exonstart = rawORF[51].split('.')
                exonstop = rawORF[52].split('.')
                startORF = int(rawORF[13])
                stopORF = int(rawORF[14])

                
                position = [int(exon) for exon in exonstart+exonstop if int(exon) >= stopORF and int(exon) <= startORF]
                position.extend([startORF, stopORF])
                position.sort(reverse=True)  # startORF, exonstop, exonstart, stopORF

                for num, pos in enumerate(position):
                    if num%2 == 0:
                        ORFstartPosition.append(pos) # 짝수 : start
                    else:
                        ORFstopPosition.append(pos) # 홀수 : stop
       
                ORFstartPosition = '.'.join(map(str, ORFstartPosition))
                ORFstopPosition ='.'.join(map(str, ORFstopPosition))
                rawORF.extend([ORFstartPosition, ORFstopPosition])
                out.write('\t'.join(map(str, rawORF)) + '\n')

