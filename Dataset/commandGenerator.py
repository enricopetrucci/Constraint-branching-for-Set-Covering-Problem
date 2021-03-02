from os import listdir
from os.path import isfile, join

new_name="run1.txt"
out = open(new_name, "w")

mypath="./ReducedLP/"#"./Rail/" #
files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
files.sort()
print(files)
for k in files:
    out.write("srun ./main.o -file ../Dataset/ReducedLP/" + k + " -callback 1 -branching 0 constraintBranchVer 1 -delta 2 -seed 0 -timelimit 3600 -threads 0 -storeResults 1\n" )
