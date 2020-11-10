import pandas as pd
from os import listdir
from os.path import isfile, join

folders = ["1_0_0_0"]#, "1_1_0_2", "1_1_1_2", "1_1_2_2"]
pdDataframes = [];

flags = [0,0,0,0]
folder = 0
for i in folders:
    mypath="./"+i
    dataframe=pd.DataFrame()
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    for f in files:
        df=pd.read_csv(mypath+"//"+f)
        dataframe = dataframe.append(df);
        
    print(dataframe)
    folder+=1
    
