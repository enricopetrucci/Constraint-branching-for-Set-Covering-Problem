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

folderDict={"0": "1_0_0_0", "1" : "1_1_0_2", "2": "1_1_1_2", "3": "1_1_2_2", "4": "1_1_3_2"}

if(options.f==str(-1)):
    print("Considering all folders")
    #folders = ["1_0_0_0", "1_1_0_2", "1_1_0_4", "1_1_1_2", "1_1_2_2", "1_1_3_2"]# "1_1_1_2", "1_1_2_2"]#[]##, "1_1_0_2", "1_1_1_2", "1_1_2_2"]
    folders = ["1_0_0_0", "1_1_3_5.0_500", "1_1_3_5_500", ]#"1_1_3_2.0_500", "1_1_3_1.0_1000", "1_1_3_1.1_500", "1_1_3_1.5_100", "1_1_3_1.5_500", "1_1_3_1.5_1000", []##, "1_1_0_2", "1_1_1_2", "1_1_2_2"] "1_1_3_2_100", "1_1_3_2.0_100", "1_1_3_2_10", "1_1_3_2_100", "1_1_3_2_500", "1_1_3_5_100", "1_1_3_2_500", "1_1_3_1.5_100","1_1_3_1.5_500", "1_1_3_1.1_500", "1_1_3_1.0_1000", "1_1_3_1.5_1000", "1_1_2_2"
else:
    folders = [folderDict[str(options.f)]] 
               
pdDataframes = []

folder = 0
for i in folders:
    mypath="./"+i
    dataframe=pd.DataFrame()
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    files.sort()
    for f in files:
        df=pd.read_csv(mypath+"//"+f)
        dataframe = dataframe.append(df,ignore_index=True);
    
    if options.seed != '-1':
        print("masking")
        mask = (dataframe['Seed'] == int(options.seed))
        dataframe = dataframe[mask]
    
    pdDataframes.append(dataframe)
    print(dataframe)
    folder+=1


maskO = (pdDataframes[0]['Time'] > int(options.timeOffset))


for i in range(len(pdDataframes)):    
    pdDataframes[i] = pdDataframes[i][maskO]    
    mask = (pdDataframes[i]['Time'] < 3600)
    pdDataframes[i] = pdDataframes[i][mask]
    print(pdDataframes[i])


##manage selecting the seed
numInstances = len(pdDataframes[0].index)

print("considering ", numInstances, " istances")

for i in pdDataframes:
    if(numInstances!=len(i.index)):
        print("Error")
        sys.exit()

##check if they are all with the same number of threads


count = 0
for i in pdDataframes:    
    timeCol = i["Time"].tolist()
    shift=int(options.timeShift)
    shiftedTimes = [x+shift for x in timeCol]

    geometricMean = gmean(shiftedTimes) - shift

    print('The geometric Mean for the runs in the ',folders[count],' folder is '+ str(geometricMean)) 
    count+=1    

    
numInstances = len(pdDataframes[0].index)

for i in pdDataframes:
    if(numInstances!=len(i.index)):
        print("Error")
        sys.exit()


count = 0
for i in pdDataframes:    
    NodesCol = i["Nodes"].tolist()
    shift=int(options.nodeShift)
    shiftedNodes = [x+shift for x in NodesCol]

    geometricMean = gmean(shiftedNodes) - shift

    print('The geometric Mean of the nodes for runs in the ',folders[count],' folder is '+ str(geometricMean)) 
    count+=1


count = 0
for i in pdDataframes:    
    branchConst = i["Constraint Branching"].tolist()
    shift=int(options.nodeShift)
    
    shiftedbranchConst = [x+shift for x in branchConst]

    geometricMeanBr = gmean(shiftedbranchConst) - shift

    print('The geometric Mean of the branching on constraints for runs in the ',folders[count],' folder is '+ str(geometricMeanBr)) 

    branchVar = i["Varaible Branching"].tolist()
    shift=int(options.nodeShift)
    
    shiftedbranchVar = [x+shift for x in branchVar]

    geometricMeanVar = gmean(shiftedbranchVar) - shift

    print('The geometric Mean of the branching on variables nodes for runs in the ',folders[count],' folder is '+ str(geometricMeanVar)) 
    count+=1
