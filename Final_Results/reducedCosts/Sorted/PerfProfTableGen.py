import pandas as pd
import sys
from os import listdir
from os.path import isfile, join
from scipy.stats.mstats import gmean 


from optparse import OptionParser

parser = OptionParser()

parser.add_option("-f","--folder", dest="f", default='-1',
                  help="select folder to consider")

parser.add_option("-r","--seed", dest="seed", default='-1',
                  help="don't print status messages to stdout")

parser.add_option("--timeShift", dest="timeShift", default='-1',
                  help="don't print status messages to stdout")

parser.add_option("--nodeShift", dest="nodeShift", default='-1',
                  help="don't print status messages to stdout")

parser.add_option("--timeOffset", dest="timeOffset", default=0,
                  help="don't print status messages to stdout")



(options, args) = parser.parse_args()

print(options)
print(options.f)

indexes = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10,]

names =["CPLEX Default", 
        "lunghezza crescente",
        "score reduced costs somma crescente",
        "score reduced costs medio crescente",
        "lunghezza crescente duplicati in cima",
        "score reduced costs somma crescente duplicati in cima",
        "score reduced costs medio crescente duplicati in cima",
        "score reduced costs somma decrescente",
        "score reduced costs medio decrescente",
        "score reduced costs somma decrescente duplicati in cima",
        "score reduced costs medio decrescente duplicati in cima"
        ]

folders = [   
            "1_0_0_0.0_0_0_0_0_0",
            "1_1_3_1.5_50_0_0_0_0",
            "1_1_3_1.5_50_0_0_2_0",
            "1_1_3_1.5_50_0_0_2_1",
            "1_1_3_1.5_50_0_1_0_0",
            "1_1_3_1.5_50_0_1_2_0",
            "1_1_3_1.5_50_0_1_2_1",
            "1_1_3_1.5_50_1_0_2_0",
            "1_1_3_1.5_50_1_0_2_1",
            "1_1_3_1.5_50_1_1_2_0",
            "1_1_3_1.5_50_1_1_2_1"
            ]

               
pdDataframes = []
notSolved = []
folder = 0

for i in indexes:
    mypath="./"+folders[i]
    dataframe=pd.DataFrame()
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files.sort()
    for f in files:
        df=pd.read_csv(mypath+"//"+f)
        dataframe = dataframe.append(df,ignore_index=True);
    
    pdDataframes.append(dataframe)
    print(dataframe)
    folder+=1

maskO = (pdDataframes[0]['Time'] > int(options.timeOffset))

for i in range(len(pdDataframes)):    
    pdDataframes[i] = pdDataframes[i][maskO]

for i in range(len(pdDataframes)):
    notSolved.append(0)
    times = pdDataframes[i]["Time"].tolist()
    for j in times:
        if(j>3600):
            notSolved[i]+=1


##manage selecting the seed
numInstances = len(pdDataframes[0].index)

for i in pdDataframes:
    if(numInstances!=len(i.index)):
        print("Error")
        sys.exit()

d = {str(len(indexes)): pdDataframes[0]['Instance']}
for i in indexes:
    d[names[i]] = pdDataframes[i]['Time']
df = pd.DataFrame(data=d)

df.to_csv("1.csv", index=False)


d = {str(len(indexes)): pdDataframes[0]['Instance']}
for i in indexes:
    d[names[i]] = pdDataframes[i]['Nodes']
df = pd.DataFrame(data=d)

df.to_csv("2.csv", index=False)
