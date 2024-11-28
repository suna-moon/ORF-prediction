### 그래프 그리기 
##: x - RibORF [PME], y - RiboCode [-log(adj P value)]  => RibORF, RiboCode의 prediction이 얼마나 예측 ORF 공통점(이라기보다는.. 유사도?...라기보다는...)을 가지는지 확인
# PME : Col47, adjusted_pval : Col30
import matplotlib.pyplot as plt
import pandas as pd
import math
min = {'D0':1.331688e-307, 'D4': 1.240419e-307, 'D8':3.349093e-294}

# data 불러오기
data= pd.DataFrame()
for i in 0, 4, 8: 
    df = pd.read_csv(f'/Data_3/Suna/ORF/commonORF/D{i}.ribocode.riborf.csv', sep='\t')
    # display(df[df['adjusted_pval'] == 0]) # 0으로 표시되는 금쪽s
    # display(df[df['adjusted_pval'] == 0].iloc[:, [0, 1, 29, 46]])

    df = df.iloc[:, [29, 46]]
    df.loc[df['adjusted_pval'] == 0, 'adjusted_pval'] = min[f'D{i}']

    display(len(df[df['adjusted_pval'] <= 3.349093e-294].iloc[:, [0, 1]])) # 금쪽이 확인용



    #### 최솟값 알아내기 목적 => 'D0':1.331688e-307, 'D4': 1.240419e-307, 'D8':3.349093e-294
    # display(df.sort_values(by='adjusted_pval', axis=0, ascending=True))

    data = pd.concat([data, df], axis=1)
display(data)
data.columns = ['D0(adjusted_pval)', 'D0(PME)', 'D4(adjusted_pval)', 'D4(PME)', 'D8(adjusted_pval)', 'D8(PME)']
display(data.sort_values(by='D0(adjusted_pval)', axis=0, ascending=True))

# 그래프
fig = plt.figure(figsize=[20, 5])
day = [0, 4, 8]

titleFsize =25
suptitleFsize = 25
label = 23

# D0(adjusted_pval)	D0(PME)	D4(adjusted_pval)	D4(PME)	D8(adjusted_pval)	D8(PME)
for i, j in enumerate(day):
    x = list(data[f'D{j}(PME)'])
    y = [math.log10(k)*(-1) for k in data[f'D{j}(adjusted_pval)']]
    ax = fig.add_subplot(1, 3, i+1)
    ax.set_xlabel('RibORF (PME)', fontdict={'fontsize':label}, labelpad=10)
    ax.set_ylabel('RiboCode\n-log(Adjusted P-value)', fontdict={'fontsize':label}, labelpad=10)
    
    # 제목
    ax.set_title(f'D{j}', fontsize=titleFsize)
    ax.scatter(x, y, s =5, c='#FF9494', alpha=0.5)


plt.suptitle('Comparison of statistical values ​​between RiboORF & RiboCode', fontsize=suptitleFsize)  
plt.tight_layout()
plt.savefig('/Data_3/Suna/ORF/commonORF/statistics(PME,adj-P).pdf', dpi=300, bbox_inches= 'tight')
plt.show()
