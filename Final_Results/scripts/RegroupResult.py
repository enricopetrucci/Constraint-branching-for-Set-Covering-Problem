import pandas as pd
import sys
from os import listdir
from os.path import isfile, join
from scipy.stats.mstats import gmean 


folders = [
            "1_0_0_0.0_0_0_0_0_0",
            "1_1_3_1.5_50_0_0_0_0",
            "1_1_3_1.5_50_0_0_1_0",
            "1_1_3_1.5_50_0_0_1_1",
            "1_1_3_1.5_50_0_1_0_0",
            "1_1_3_1.5_50_0_1_1_0",
            "1_1_3_1.5_50_0_1_1_1",
            "1_1_3_1.5_50_1_0_0_0",
            "1_1_3_1.5_50_1_0_1_0",
            "1_1_3_1.5_50_1_0_1_1",
            "1_1_3_1.5_50_1_1_0_0",
            "1_1_3_1.5_50_1_1_1_0",
            "1_1_3_1.5_50_1_1_1_1"
            ]

pdDataframes = []
folder = 0
for i in folders:
    dataframe=pd.DataFrame()
    files = [f for f in listdir("./Final_result5/"+i) if isfile(join("./Final_result5/"+i, f))]
    files.sort()
    for f in files:
        df0 = pd.read_csv("./Final_result5/"+i+"//"+f)
        df1 = pd.read_csv("./Final_result5Alternative/"+i+"//"+f)
        #df0 = df0.append(df1,ignore_index=True);
        df2 = pd.read_csv("./Final_result5Alternative1/"+i+"//"+f)
        df0 = df1.append(df2,ignore_index=True);
              
        df0.to_csv("./10Seed/"+i+"/"+f, index=False)
