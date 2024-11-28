import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np


# 함수
def threediv(graphdf, stop):
    data_0 = pd.DataFrame(columns=[0, 1])
    data_1 = pd.DataFrame(columns=[0, 1])
    data_2 = pd.DataFrame(columns=[0, 1])
    
    for i in graphdf.index:
        if i % 3 == 0:
            data_0.loc[len(data_0)] = [i, graphdf.at[i, stop]]
        elif i%3 == 1:
            data_1.loc[len(data_1)] = [i, graphdf.at[i, stop]]
        else:
            data_2.loc[len(data_2)] = [i, graphdf.at[i, stop]]
    return data_0, data_1, data_2

def threent (df, start):
    zero = 0
    one = 0
    two = 0
    for i in df[start]:
        if i % 3 == 0:
            zero += 1
        elif i % 3 == 1 :
            one += 1
        else: # i % 3 == 2 :
            two += 1
    dfCount = pd.DataFrame({0:[zero], 1:[one], 2:[two]})
    return dfCount

# 전체 gene, 3개 그래프
top_list = ['Tshr', 'Dgat2', 'Ppp2r5a', 'Ucp2', 'Hr', 'Acp6', 'Mrpl42', 'Abcb6', 'Crls1', 'Ebf2', 'Btg1', 'Slc25a44',
  'Ppp2r5b', 'Ppp1r15a', 'Ptpa', 'Stambpl1', 'Ppp1r15b', 'Gnpat', 'Timm13', 'Pcnx4', 'Mief1', 'Tor1aip1', 'Ptp4a1_0', 'Ptp4a1_1', 'Kdm5c', 'Klhl24']
bottom_list = ['Slit2', 'Fibin', 'Cdh11', 'Prss23', 'Vgll3', 'Igfbp4', 'Col3a1', 'Slit3', 
               'Pcdhgc3', 'Rhoj', 'Mgat2', 'Pkd2', 'Gas6', 'Lrrc8a', 'Pik3r2', 'Maf', 'Mbtps1', 'Mlec', 'Loxl1']
total_list = top_list+bottom_list
for gene_name in total_list:
    for x in [0, 4, 8]:
        for y in ['a', 'b', 'c']:
            rep = f'D{x}{y}'
            filepath = f'/Data_3/Suna/ORF/3nt/count/{gene_name}_{rep}.txt'
            df = pd.read_csv(filepath, sep='\t', header=0)
            start = df.columns[0]
            stop = df.columns[1]
            if x == 0 and y == 'a':  # 처음 loop일 때 # down 부분 아래 1개로 합치기
                myArr = np.zeros((int(stop)-int(start)+1, 1)) # start 부터 stop까지 0인 sum df지정
                sum = pd.DataFrame(myArr, columns=[stop], index=[str(i) for i in np.arange(int(start), int(stop)+1, 1)])
            dfforsum = df.copy() # index를 변경
            dfforsum[start] = dfforsum[start].astype(str)
            dfforsum = dfforsum.set_index(start)
            sum = dfforsum.add(sum, fill_value=0)

        df = sum.copy(deep=True).reset_index()
        df.columns = [start, stop]
        df[start] = df[start].astype(int)-int(start) # start 0 만들어주기


        # norm
        norm = []
        k = int(max(df[start]))
        normdf = pd.DataFrame()
        graphdf = pd.DataFrame()
        df = df.set_index(start)

        # read 개수 매우 적을 때 pd.concat 안되는 경우
        if k % 3  == 0:
            k = k + 1 + 2
        elif k % 3 == 1 :
            k = k + 1 + 1
        else:
            k = k + 1 + 0

        # start~stop norm
        for i in range(0, k):
            if i % 3 == 2:
                try:
                    normdf[i] = df.loc[i]
                except KeyError: # str '-'일 때 0
                    normdf.at[stop, i] = 0
                
                test = normdf.sum(axis=1).values
                # print(test, type(test[0]))
                if normdf.sum(axis=1).values == 0: # 0, 1, 2 다 돌고 합이 9이면 3개 0, 0, 0
                    normdf.loc[stop] = [0, 0, 0]
                else:
                    normdf = normdf/normdf.sum(axis=1).values[0]*100 # 전체 read 수로 norm

                graphdf = pd.concat([graphdf, normdf], axis=1) # 값전달
                normdf=pd.DataFrame() #초기화
            else:
                try:
                    normdf[i] = df.loc[i]
                except KeyError:
                    normdf.at[stop, i] = 0

        ## plot
        if x == 0:
            graphdf_D0 = graphdf.T.copy(deep=True)
        elif x== 4:
            graphdf_D4 = graphdf.T.copy(deep=True)
        else:
            graphdf_D8 = graphdf.T.copy(deep=True)
    
    # fig = plt.figure(figsize=(30,12)) # 36, 12
    fig = plt.figure(figsize=(10,3)) # 36, 12
    # gs = gridspec.GridSpec(nrows=3, ncols=2, width_ratios=[10, 2], height_ratios=[1, 1, 1])
    # gs = gridspec.GridSpec(nrows=3, ncols=1, height_ratios=[1, 1, 1])
    gs = gridspec.GridSpec(nrows=1, ncols=1) # 1개

    ax5 = plt.subplot(gs[0])

    # 그래프 data
    D00, D01, D02 = threediv(graphdf_D0, stop)
    D40, D41, D42 = threediv(graphdf_D4, stop)
    D80, D81, D82 = threediv(graphdf_D8, stop)

    wid = 0.8
    linewid = 1

    ax5.bar(D80[0], D80[1], width=wid, zorder=1, linewidth=linewid, edgecolor='white', color=['lightcoral'])
    ax5.bar(D81[0], D81[1], width=wid, zorder=1, linewidth=linewid, edgecolor='white', color=['silver'])
    ax5.bar(D82[0], D82[1], width=wid, zorder=1, linewidth=linewid, edgecolor='white', color=['grey'])

    font = 10

    ax5.set_ylabel('Frame Percent\nWithin Codons(%)', fontsize=font)
    ax5.set_xticks([0, int(stop)-int(start)+1])
    ax5.set_xticklabels(['uORF\nStart Codon', 'uORF\nStop Codon'], fontsize=font)

    ax5.set_ylim(0, 100)

    b = 1

    ax5.set_xlim(0-b, int(stop)-int(start) + 3.5)

    ax5.spines['top'].set_visible(False) #오른쪽 테두리 제거
    ax5.spines['right'].set_visible(False) #오른쪽 테두리 제거

    # 제목
    ax5.set_title(f'{gene_name}', pad=20, fontsize=font)

    # # 범례
    ax5.legend(['Frame 0', 'Frame 1', 'Frame 2'], loc = 'lower center', bbox_to_anchor = (0,-0.3,1,1), ncols = 3, edgecolor='white', fontsize=font)

    # Show graph
    print(gene_name, start, stop, int(stop)-int(start))
    plt.tight_layout()
    plt.savefig(f"/Data_3/Suna/ORF/3nt/fig/one/{gene_name}.pdf", dpi = 300, bbox_inches = "tight")
    plt.show()
