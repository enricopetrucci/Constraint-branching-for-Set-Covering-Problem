import pandas as pd 
import sys 
from os import listdir 
from os.path import isfile, join 
from scipy.stats.mstats import gmean  
 
import matplotlib 
import matplotlib.pyplot as plt 
import numpy as np 
 
folders = [ 
        "1_2_0_0.0_0_0_0_0_0", 
        "1_3_0_1.0_0_0_0_0_0", 
        "1_3_1_1.0_0_0_0_0_0" 
        ] 
 
 
instances = [ 
"scpb2_reduced", 
"scpb4_reduced", 
 
"scpc2_reduced", 
"scpc3_reduced", 
 
"scpd1_reduced", 
"scpd2_reduced", 
"scpd3_reduced", 
"scpd4_reduced", 
 
#"scpnre1_reduced", 
#"scpnre2_reduced", 
#"scpnre3_reduced", 
 
#"scpnrf1_reduced", 
#"scpnrf2_reduced", 
#"scpnrf3_reduced", 
] 
 
seeds = [338, 1867, 6822, 8927, 9397] 
 
files = [] 
for i in instances: 
    for j in seeds: 
        files.append(i + "_" + str(j) + ".csv") 
         
print(files) 
print("considering ", len(files), " computations") 
 
g = [] 
var = [] 
t = [] 
 
for k in folders: 
    mypath="./"+k 
    print("considering folder = ", k) 
    dataframe=pd.DataFrame() 
    for f in files: 
        df=pd.read_csv(mypath+"/"+f) 
        dataframe = dataframe.append(df,ignore_index=True); 
         
    print(dataframe) 
 
    print("maximum depth = ", dataframe.max(axis=0)["depth"]) 
 
    averageFracVarPerc = [] 
    averageGap = [] 
    for i in range(dataframe.max(axis=0)["depth"]): 
        mask = dataframe['depth'] == i 
        dfCurrentDepth = dataframe[mask] 
        print(dfCurrentDepth) 
        variables = dfCurrentDepth["percentageFracVars"].tolist() 
        shift = 0.01 
        shiftedVariables = [x+shift for x in variables] 
        averageFracVarPerc.append(gmean(shiftedVariables) - shift) 
 
        gaps = dfCurrentDepth["Gap"].tolist() 
        shift = 0.01 
        shiftedgaps = [x+shift for x in gaps] 
        averageGap.append(gmean(shiftedgaps) - shift) 
 
    t.append(np.arange(0, dataframe.max(axis=0)["depth"], 1)) 
 
    var.append(averageFracVarPerc) 
    g.append(averageGap) 
    print(averageFracVarPerc) 
    print(averageGap) 
 
 
 
# Data for plotting 
 
s0 = var[0] 
s1 = var[1] 
s2 = var[2] 
 
fig, ax = plt.subplots() 
ax.plot(t[0], s0, label="CPLEX Full Strong Branching") 
ax.plot(t[1], s1, label="Full Strong Constraint Branching") 
ax.plot(t[2], s2, label="Strong Constraint Branching") 
plt.legend() 
 
 
ax.set(xlabel='Depth', ylabel='Percentage of fractional variables', 
       title='Average percentage of frac. variables w.r.t branching tree depth') 
ax.grid() 
 
fig.savefig("percentageFracVars.pdf") 
plt.show() 
 
 
#Data for plotting 
s0 = g[0] 
s1 = g[1] 
s2 = g[2] 
 
fig, ax = plt.subplots() 
ax.plot(t[0], s0, label="CPLEX Full Strong Branching") 
ax.plot(t[1], s1, label="Full Strong Constraint Branching") 
ax.plot(t[2], s2, label="Strong Constraint Branching")  
plt.legend() 
 
ax.set(xlabel='Depth', ylabel='Gap', 
       title='Average Gap w.r.t branching tree depth') 
ax.grid() 
 
fig.savefig("Gap.pdf") 
plt.show() 
