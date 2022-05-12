import pandas as pd
import numpy as np
import math

info = pd.read_csv('info.csv')
# print(info)
SH = info[["上海"]]
AH = info[["安徽"]]
JS = info[["江苏"]]
FJ = info[["福建"]]
ZJ = info[["浙江"]]
# print(ZJ)
SH_new = np.zeros(math.ceil(SH.shape[0]/4))
AH_new = np.zeros(math.ceil(SH.shape[0]/4))
JS_new = np.zeros(math.ceil(SH.shape[0]/4))
FJ_new = np.zeros(math.ceil(SH.shape[0]/4))
ZJ_new = np.zeros(math.ceil(SH.shape[0]/4))
# date = np.zeros(math.ceil(SH.shape[0]/4))
# time = np.zeros(math.ceil(SH.shape[0]/4))
date = []
time = []
print(SH.shape[0])
print(info.at[0,"上海"])
for i in range(SH.shape[0]):
    if i%4 == 0:
        SH_new[int(i/4)] = (info.at[i,"上海"]+info.at[i+1,"上海"]+info.at[i+2,"上海"]+info.at[i+3,"上海"])/4
        AH_new[int(i/4)] = (info.at[i,"安徽"]+info.at[i+1,"安徽"]+info.at[i+2,"安徽"]+info.at[i+3,"安徽"])/4
        JS_new[int(i/4)] = (info.at[i,"江苏"]+info.at[i+1,"江苏"]+info.at[i+2,"江苏"]+info.at[i+3,"江苏"])/4
        FJ_new[int(i/4)] = (info.at[i,"福建"]+info.at[i+1,"福建"]+info.at[i+2,"福建"]+info.at[i+3,"福建"])/4
        ZJ_new[int(i/4)] = (info.at[i,"浙江"]+info.at[i+1,"浙江"]+info.at[i+2,"浙江"]+info.at[i+3,"浙江"])/4
        date.append(info.at[i,"date"])
        time.append(info.at[i,"time"])
print(SH_new)

dict = {"date":date,'time':time,"ShangHai":SH_new,"AnHui":AH_new,"JiangSu":JS_new,"FuJian":FJ_new,"ZheJiang":ZJ_new }
df = pd.DataFrame(dict)
df.to_csv('reslut.csv')