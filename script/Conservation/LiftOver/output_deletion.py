import pandas as pd


def skip_logic(index, skip_num):
    if index % skip_num == 1:
        return False
    return True

cate = ['uORF', 'dORF', 'lncRNA']
for i in cate:
    input_df = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/{i}.bed6', sep='\t', header=None)
    del_df = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/output/deletion/{i}_deletion.bed', sep='\t', header=None, skiprows=lambda x: skip_logic(x, 2))
    
    ### deletion 제거(예측 안된 것)
    merg_df = pd.concat([input_df, del_df])
    merg_df.drop_duplicates(keep=False, inplace=True)

    ### mouse + human
    human = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/output/{i}_liftover.bed', sep='\t', header=None)
    mouse_human = pd.concat([merg_df.reset_index(drop=True), human], axis=1, ignore_index=True)
    mouse_human.to_csv(f'{i}_mouse_human.txt', sep='\t', header=None, index=None)
display(mouse_human)

### human 파씽 : /Data_3/Suna/parsing/parsed_human.txt
cate = ['uORF', 'dORF', 'lncRNA']
# for i in cate:
ref_df = pd.read_csv('/Data_3/Suna/parsing/parsed_human.txt', sep='\t')
ref_df = ref_df.iloc[:, [0, 3, 4, 5, 11, 12, 13, 14]]  # gene_name, strand, seqname(chr?), exon/CDS+start/stop
exon_cordinate = []
write_cordinate = []
start = ref_df.exon_start.str.split('.')
stop = ref_df.exon_stop.str.split('.')
for i, j in zip(start, stop):
        for k, l in zip(i, j):
                exon_cordinate.append([k, l])
        write_cordinate.append(exon_cordinate)
        exon_cordinate = []  
ref_df['exon_cordinate'] = write_cordinate
display(ref_df)
## repre만 남기기, #/Data_1/References/hg38_v42/representative-isoforms.txt 
repre_df = pd.read_csv('/Data_1/References/hg38_v42/representative-isoforms.txt', sep='\t', header=None)

ref_df = pd.merge(ref_df, repre_df, how='inner',
                left_on='transcript_id', right_on=1).drop(labels=['transcript_id', 0, 1, 2, 3, 4], axis=1)
display(ref_df)    


cate = ['uORF', 'dORF', 'lncRNA']
for i in cate:
    with open(f'/Data_3/Suna/ORF/Liftover/compare/{i}_mouse_human.txt', 'r') as f, open(f'/Data_3/Suna/ORF/Liftover/compare/human_cordinate/{i}_human_cordinate.txt', 'w') as out:
        for line in f:
            data = line.rstrip('\n').split('\t')
            for chr, group in ref_df.groupby('seqname'):
                if chr == data[6]:  # parsing chr와 liftover data의 chr 같은 것에 대하여
                    # 이제 exon 봐야지
                    for cordinate, name in zip(group['exon_cordinate'], group['gene_name']):
                        info = []
                        info.append(name)
                        info.append(chr)
                        for exon in cordinate:
                            if exon[0] <= data[7] and data[8] <= exon[1]:
                                if name == data[9].upper():  ####### gene_name 비교
                                    write_info = []
                                    write_info.extend(data) # liftover 정보
                                    write_info.extend(info) # parsing gene_name, chromosome 정보
                                    write_info.extend(exon) # parsing 그중 exon 정보
                                    # data = liftover, parsing의 data를 가지고 옴
                                    out.write('\t'.join(write_info) + '\n') ### data 입력
                                    break
                else:
                    continue
                  
### 중복제거
cate = ['uORF', 'dORF', 'lncRNA']
for i in cate:
    df = pd.read_csv(f'/Data_3/Suna/ORF/Liftover/compare/human_cordinate/{i}_human_cordinate.txt', sep='\t', header=None)
    df.drop_duplicates(subset=[12, 14, 15], keep=False, inplace=True)
    df.to_csv(f'/Data_3/Suna/ORF/Liftover/compare/dup/{i}_human_cordinate_duplicate.txt', sep='\t', header=None, index=None)
