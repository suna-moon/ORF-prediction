import matplotlib.pyplot as plt
import pandas as pd
countDF = pd.DataFrame()
df2 = pd.DataFrame()
j = 0
for i in 0, 4, 8: 
    
    df = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{i}.ribocode.riborf.csv', sep='\t')
    df = df.groupby('candidateORFType').count().iloc[:, [0]].T # data ORF type 별로 count
    df['sum'] = df.sum(axis=1) # sum 계산
    countDF = pd.concat([countDF, df])
    countDF.index.values[j] = f'D{i}'
    j +=1
display(countDF)



#### df2 countDF3에서 transcript_type 분석하기 위함
for i in 0, 4, 8: 
    
    df = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{i}.ribocode.riborf.csv', sep='\t')
    day = f'D{i}'
    df.insert(51, 'day', [day]*len(df), True)    
    df2 = pd.concat([df2, df])

df2=df2[df2['candidateORFType'] == 'noncoding']
# .groupby(['day', 'transcript_type']).count().iloc[:,[0]]
# test file
# df2 = df2.drop_duplicates(['ORF_type', 'transcript_type', 'candidateORFType'], keep='first').loc[:, ['ORF_type', 'transcript_type', 'candidateORFType']]
display(df2)
# percentage
countDF1 = countDF.copy(deep=True) 
for i in range(3):
    countDF1.iloc[i] = countDF.iloc[i]/countDF.values[i][9]*100

# canonical vs noncanonical ORF % 계산
countDF5 = pd.DataFrame()
countDF5['canonical'] = countDF1.loc[:, 'canonical']
countDF5['noncanonical'] = 100-countDF1.loc[:, 'canonical']
countDF5 = countDF5.sort_values(by=['D0', 'D4', 'D8'], axis=1, ascending=False)


######## 2번째 non-canonical 중에서 uORF, dORF, ncORF(noncoding), others(dORF, iORF, odORF, ouORF, truncation, extension)
countDF2 = countDF.drop(['canonical', 'sum'], axis=1)
countDF2['sum'] = countDF2.sum(axis=1)
display(countDF2)
for i in range(3):
    countDF2.iloc[i] = countDF2.iloc[i]/countDF2.values[i][8]*100
countDF2 = countDF2.drop(['sum'], axis=1)
countDF2 = countDF2.sort_values(by=['D0', 'D4', 'D8'], axis=1, ascending=False)
countDF2.rename(columns={'noncoding':'ncORF'}, inplace=True)
display(countDF2)


# ######## 3번째 noncoding에 대한 종류m## 3번째 graph ncORF의 type(trnascrip_type)별 분포## 3번째 graph data : ncORF의 type(trnascrip_type)별 분포
countDF3 = pd.DataFrame()
for i in 0, 4, 8:
    countDF3 = pd.concat([countDF3, df2[df2['day']==f'D{i}'].groupby('transcript_type').count().iloc[:, [0]].T]).rename(index={'ORF_ID':f'D{i}'})
countDF3.fillna(value=0, inplace=True)
display(countDF3)
day = [0, 4, 8]
df3 = countDF3.copy(deep=True)
countDF3 = pd.DataFrame([countDF3.loc[f'D{j}']/countDF3.sum(axis=1)[i]*100 for i, j in enumerate(day)]).sort_values(by=['D0', 'D4', 'D8'], axis=1, ascending=False)
display(countDF3)

# 그래프 그리기

# fontsize
suptitle = 23
title = 15
xlabelFsize = 10
ylabelFsize = 15
legend = 15
padsize = 10
textFsize = 10



### ax1에 %값 추가 ax2에 i, overlap 애들도 추가

fig = plt.figure(figsize=[12, 5])
ax1 = fig.add_subplot(1, 2, 1)
ax2 = fig.add_subplot(1, 2, 2)
# ax3 = fig.add_subplot(1, 3, 3)

####### 1번째 canonical vs noncanonical ORF관련 % 계산
countDF5.plot.bar(stacked=True, rot=0, cmap='Pastel2_r', ax=ax1) ### 색조정필요
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles[::-1], labels[::-1], loc='lower left', bbox_to_anchor=(1, 0.5), frameon=False, fontsize=15)
ax1.set_title('Canonical vs Noncanonical', fontsize=title, pad=padsize)
ax1.set_ylim([0, 100])
ax1.set_ylabel('Percentage (%)', fontdict={'fontsize':ylabelFsize})
ax1.set_xlabel('Adipocytes Differentiation (Day)', fontdict={'fontsize':xlabelFsize})
# text % 표시해주기

display(countDF5)
values1 = list(countDF5['canonical'])
values2 = list(countDF5['noncanonical'])
for i, values1 in enumerate(values1) : 
  ax1.text(i, values1+values2[i]/2-1.5, f'{round(values2[i])}%', ha = 'center', weight = 'bold', color = 'black', fontdict={'fontsize':textFsize})
  ax1.text(i, values1/2, f'{round(values1)}%', ha = 'center', weight = 'bold', color = 'black', fontdict={'fontsize':textFsize})

######## 2번째 non-canonical 중에서 uORF, dORF, ncORF(noncoding), others(dORF, iORF, odORF, ouORF, truncation)
countDF2.plot.bar(stacked=True, rot=0, cmap='Accent', ax=ax2)
handles, labels = ax2.get_legend_handles_labels()
display(countDF2)
# a = list(reversed(list(countDF2.columns)))
ax2.legend(handles[::-1], labels[::-1], loc='lower left', bbox_to_anchor=(1, 0.1), frameon=False, fontsize=15)
ax2.set_title('Noncanonical ORF Classification', fontsize=title, pad=padsize)
ax2.set_ylim([0, 100])
ax2.set_ylabel('Percentage (%)', fontdict={'fontsize':ylabelFsize})
ax2.set_xlabel('Adipocytes Differentiation (Day)', fontdict={'fontsize':xlabelFsize})

######## 3번째 noncoding에 대한 종류m## 3번째 graph ncORF의 type(trnascrip_type)별 분포
# countDF3.plot.bar(stacked=True, rot=0, cmap='Pastel2', edgecolor='white', ax=ax3)
# handles, labels = ax3.get_legend_handles_labels()
# ax3.legend(handles[::-1], labels[::-1], loc='lower left', bbox_to_anchor=(1, 0.25), frameon=False)
# ax3.set_title('Transcript Type of ncORF', fontsize=title, pad=padsize)
# ax3.set_ylim([0, 100])
# ax3.set_ylabel('Percentage (%)', fontdict={'fontsize':labelFsize})
# ax3.set_xlabel('Adipocytes Differentiation (Day)', fontdict={'fontsize':labelFsize})


fig.suptitle('ORF Ratio', fontsize=suptitle)
plt.tight_layout()
plt.savefig(f"/Data_3/Suna/ORF/commonORF/ORFratio.pdf", dpi = 300, bbox_inches = "tight")
plt.show() 


