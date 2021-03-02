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
    folders = [
               "1_0_0_0.0_0_0_0_0_0",
               
               "1_1_3_1.5_50_0_1_0_0",
               "1_1_3_1.5_100_0_1_0_0",              
               "1_1_3_1.5_150_0_1_0_0",
               "1_1_3_1.5_200_0_1_0_0",
               
               "1_1_3_2.0_50_0_1_0_0",
               "1_1_3_2.0_100_0_1_0_0",
               "1_1_3_2.0_150_0_1_0_0",
               "1_1_3_2.0_200_0_1_0_0",
               
               "1_1_3_5.0_50_0_1_0_0",
               "1_1_3_5.0_100_0_1_0_0",
               "1_1_3_5.0_150_0_1_0_0",
               "1_1_3_5.0_200_0_1_0_0",
               
               ]#"1_1_3_2.0_500", "1_1_3_1.0_1000", "1_1_3_1.1_500", "1_1_3_1.5_100", "1_1_3_1.5_500", "1_1_3_1.5_1000", []##, "1_1_0_2", "1_1_1_2", "1_1_2_2"] "1_1_3_2_100", "1_1_3_2.0_100", "1_1_3_2_10", "1_1_3_2_100", "1_1_3_2_500", "1_1_3_5_100", "1_1_3_2_500", "1_1_3_1.5_100","1_1_3_1.5_500", "1_1_3_1.1_500", "1_1_3_1.0_1000", "1_1_3_1.5_1000", "1_1_2_2"
else:
    folders = [folderDict[str(options.f)]] 
               
pdDataframes = []
notSolved = []
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
        dataframe.reset_index(inplace=True)
    
    #mask1 = dataframe["Instance"] != "../Dataset/ReducedLP/scpclr11_reduced.lp"
    #mask2 = dataframe["Instance"] != "../Dataset/ReducedLP/scpnrg2_reduced.lp"
    #mask = [a and b for a, b in zip(mask1, mask2)]
    #dataframe = dataframe[mask]
    #dataframe.reset_index(inplace=True)
    
    
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

print("considering ", numInstances, " istances")

for i in pdDataframes:
    print(len(i.index))
    if(numInstances!=len(i.index)):
        print("Error")
        sys.exit()

##check if they are all with the same number of threads


time = []
nodes = []
consbr = []
varbr = []

count = 0
for i in pdDataframes:    
    timeCol = i["Time"].tolist()
    shift=int(options.timeShift)
    shiftedTimes = [x+shift for x in timeCol]

    geometricMean = gmean(shiftedTimes) - shift
    time.append(geometricMean)
    print('The geometric Mean for the runs in the ',folders[count],' folder is '+ str(geometricMean)) 
    count+=1    

    
numInstances = len(pdDataframes[0].index)

for i in pdDataframes:
    print(len(i.index))
    if(numInstances!=len(i.index)):
        print("Error")
        sys.exit()


count = 0

for i in pdDataframes:
    mask = (i['Time'] < 3600)
    print("considering dataframe ", count)
    deleating = 0
    for k in mask:
        if(k == False):
            deleating+=1
    print("deleating ", deleating, " instances")
    count+=1
    if deleating > 0:
        for j in range(len(pdDataframes)):
            print("Dataframe ", j," from ", len(pdDataframes[j].index))
            pdDataframes[j] = pdDataframes[j][mask]
            #pdDataframes[j].reset_index(inplace=True)
            print("Dataframe ", j," to ",len(pdDataframes[j].index))
        

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
    nodes.append(geometricMean)

    print('The geometric Mean of the nodes for runs in the ',folders[count],' folder is '+ str(geometricMean)) 
    count+=1


count = 0
for i in pdDataframes:    
    branchConst = i["Constraint Branching"].tolist()
    shift=int(options.nodeShift)
    
    shiftedbranchConst = [x+shift for x in branchConst]

    geometricMeanBr = gmean(shiftedbranchConst) - shift
    consbr.append(geometricMeanBr)

    print('The geometric Mean of the branching on constraints for runs in the ',folders[count],' folder is '+ str(geometricMeanBr)) 

    branchVar = i["Varaible Branching"].tolist()
    shift=int(options.nodeShift)
    
    shiftedbranchVar = [x+shift for x in branchVar]

    geometricMeanVar = gmean(shiftedbranchVar) - shift
    varbr.append(geometricMeanVar)

    print('The geometric Mean of the branching on variables nodes for runs in the ',folders[count],' folder is '+ str(geometricMeanVar)) 
    count+=1

count = 0
for i in notSolved:
    print('Number of non solved instances in the ',folders[count],' folder is '+ str(i)) 
    count+=1




fields=["Method", "Time", "Nodes", "Var. B.", "Constr. B.", "Not Solved"]

for j in range (len(fields)):
    if(j<len(fields)-1):
        print(fields[j], end =" & ")
    else:
        print(fields[j], end =" \\\\\n")
        
for i in range(len(folders)):
    for j in range (len(fields)):
        if j==0:
            print(folders[i], end =" & ")
            
        elif j==1:
            print("{:10.3f}".format(time[i]), end =" & ")
            

        elif j==2:
            print("{:10.1f}".format(nodes[i]), end =" & ")
            
        elif j==3:
            print("{:10.1f}".format(varbr[i]), end =" & ")
            
        elif j==4:
            print("{:10.1f}".format(consbr[i]), end =" & ")
            
        elif j==5:
            print("{:10.1f}".format(notSolved[i]), end =" \\\\\n")
                
    
for i in range(len(pdDataframes)):
    print(folders[i] + "," + str(time[i]) + "," + str(nodes[i]) + "," + str(consbr[i]) + "," + str(varbr[i]) + "," + str(notSolved[i]))
