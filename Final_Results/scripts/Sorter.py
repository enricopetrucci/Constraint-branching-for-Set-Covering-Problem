import pandas as pd
import sys
from os import listdir
from os.path import isfile, join
from scipy.stats.mstats import gmean 


folders = [
            "1_0_0_0.0_0_0_0_0_0",
            "1_1_3_5.0_100_0_0_0_0",
            "1_1_3_5.0_100_0_0_2_0",
            "1_1_3_5.0_100_0_0_2_1",
            
            "1_1_3_5.0_100_0_1_0_0",
            "1_1_3_5.0_100_0_1_2_0",
            "1_1_3_5.0_100_0_1_2_1",
            "1_1_3_5.0_100_1_0_2_0",
            
            "1_1_3_5.0_100_1_0_2_1",
            "1_1_3_5.0_100_1_1_2_0",
            "1_1_3_5.0_100_1_1_2_1",
            
            
            ]

pdDataframes = []
notSolved = []
folder = 0
for i in folders:
    print("folder ", i)
    mypath="./"+i
    dataframe=pd.DataFrame()
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files.sort()
    for f in files:
        df=pd.read_csv(mypath+"//"+f)
        df.sort_values(by=['Seed'], inplace=True)
        df.to_csv("./Sorted/"+i+"/"+f, index=False)
