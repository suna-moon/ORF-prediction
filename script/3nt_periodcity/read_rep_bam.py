import pysam
import pandas as pd

# correction table
correction = pd.DataFrame({27:['-',12,12,'-','-','-','-','-','-',],
                           28:[12, 12, 12, 12, 12, 12, 12, 12, 12], 
                           29:[12, 12, 12, 12, 12, 12, 12, 12, 12],
                           30:[12, 12, 12, 12, 12, 12, 12, 12, 12],
                           31:[13, 13, 13, 13, 13, 13, 13, 13, 13],
                           32:[13, 13, 13, 13, 13, 13, 12, 12, 13],
                           33:[13, 13, 13, 13, 13, 13, 13, 13, 13],
                           34:[13, 13, 13, 13, 13, 13, 13, 13, 13]}, index=['D0a', 'D0b', 'D0c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c'])

# correction = pd.DataFrame({27:['-',12,12,'-','-','-','-','-','-',],
#                            28:[12, 12, 12, 12, 12, 12, 12, 12, 12], 
#                            29:[12, 12, 12, 12, 12, 12, 12, 12, 12],
#                            30:[12, 12, 12, 12, 12, 12, 12, 12, 12],
#                            31:[13, 13, 13, 13, 13, 13, 13, 13, 13],
#                            32:[13, 13, 13, 13, 13, 13, 13, 13, 13],
#                            33:[13, 13, 13, 13, 13, 13, 13, 13, 13],
#                            34:[13, 13, 13, 13, 13, 13, 13, 13, 13]}, index=['D0a', 'D0b', 'D0c', 'D4a', 'D4b', 'D4c', 'D8a', 'D8b', 'D8c'])


display(correction)
# name =[]
# for i in [0, 4, 8]:
#     for j in ['a', 'b', 'c']:
#         name.append(f'D{i}{j}')
# print(name)

#     27, 28, 29, 30, 31, 32, 33, 34
# D0a  -, 12, 12, 12, 13, 13, 13, 13
# D0b 12, 12, 12, 12, 13, 13, 13
# D0c 12, 12, 12, 12, 13, 13, 13
# D4a  -, 12, 12, 12, 13, 13, 13
# D4b  -, 12, 12, 12, 13, 13, 13
# D4c  -, 12, 12, 12, 13, 13, 13
# D8a  -, 12, 12, 12, 13, 12, 13
# D8b  -, 12, 12, 12, 13, 12, 13
# D8c  -, 12, 12, 12, 13, 13, 13

def count_elements(lst):
    counts = {}
    for item in lst:
        if item in counts:
            counts[item] += 1
        else:
            counts[item] = 1
    return counts

for i in [0, 4, 8]:
    for j in ['a', 'b', 'c']:
        file_name = f'D{i}{j}'
        print(file_name)
    # file_name = 'D0a'
        filepath = f'/Data_2/Jun/Adipocytes/tr-aln/novaseq/{file_name}.rep.bam'
        bamfile = pysam.AlignmentFile(filepath, 'rb')

        orf_list = pd.read_csv('/Data_3/Suna/ORF/3nt/transcript_convert.txt', sep='\t', header=0)
        tempwrite = {}
        ptp4a1 = 0
        colname = {}

        for row in orf_list.itertuples():
            id = row[1]
            if row[11]-20 <=1:
                start = 1
            else:
                start = row[11] - 20
            stop = row[12]
            pos = []
            for read in bamfile.fetch(id, start, stop):

                # length correction
                length = read.reference_length
                try:
                    cor_val= correction.at[file_name,length]
                except KeyError:
                    continue
                if type(cor_val) == str:
                    continue

                if read.reference_start >=start and read.reference_start <= stop:
                    pos.append(read.reference_start + 1 + cor_val)  # length 정보
                # tempwrite[row[3]].append(read.reference_start())

            # Ptp4a1 특별취급
            if row[3] == 'Ptp4a1':
                tempwrite[row[3]+f'_{ptp4a1}'] = pos
                colname[row[3]+f'_{ptp4a1}'] = [row[11], row[12]]
                ptp4a1 = ptp4a1 +1
            else:
                tempwrite[row[3]] = pos
                colname[row[3]] = [row[11], row[12]]

        print(tempwrite)
        print(colname)
        for key, val in tempwrite.items():
            column_name = colname[key] # start stop 정보
            df = pd.DataFrame.from_dict(count_elements(val), orient='index').reset_index()
            try:
                df.columns = column_name
            except ValueError: # count가 0일 때
                df = pd.DataFrame({0:[0], 1:[0]})
                df.columns = column_name
            # df.columns = column_name
            df.to_csv(f'/Data_3/Suna/ORF/3nt/count/{key}_{file_name}.txt', sep='\t', index=False)
        # print(len(tempwrite))
        bamfile.close()

        # print(list(orf_list['gene_name']))
        # ['Col3a1', 'Becn1', 'Ccnk', 'Marchf3', 'Fam214b', 'Pkd2', 'Tshr', 'Adamts12', 'Ucp2', 'Prss23', 'Loxl1', 'Maf', 'Mlec', 'Ppp1r15b', 'Btg1', 'Slc35a4', 'Gnpat', 'Abca1', 'Slit3', 'Pcnx4', 'Ptp4a1', 'Ptp4a1', 'Vgll3', 'Slc25a44', 'Lrrc8a', 'Gm28035', 'Timm13', 'Mgat2', 'Klhl24', 'Hdac5', 'Pik3r2', 'Tor1aip1', 'Mrpl42', 'Usp16', 'Bmpr2', 'Dnajc2', 'Cdh11', 'Rhoj', 'Yif1a', 'Pip5k1a', 'Snrpn', 'Fibin', 'Pcdhgc3', 'Ppp1r15a', 'Snx13', 'Ppp2r5b', 'Arl5a', 'S1pr3', 'S1pr3', 'Slit2', 'Cdipt', 'Crls1', 'Ebf2', 'Adamts1', 'Ptpa', 'Dynll1', 'Brd2', 'Mief1', 'Igfbp4', 'Trip11', 'Hr', 'Phkb', 'Cggbp1', 'Mbtps1', 'Sf3b3', 'Gtf2h1', 'Dgat2', 'Mbtps2', 'Stt3a', 'Ppp2r5a', 'Abcb6', 'Stambpl1', 'Herpud2', 'Ehbp1l1', 'Gas6', 'Acp6', 'Azi2', 'Golph3', 'Trim41', 'Tmem168', 'Kdm5c', 'Aff4']







