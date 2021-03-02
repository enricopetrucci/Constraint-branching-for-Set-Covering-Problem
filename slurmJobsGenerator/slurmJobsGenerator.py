 
files=[

"../Dataset/ReducedLP/scpb2_reduced.lp",
"../Dataset/ReducedLP/scpb4_reduced.lp",

"../Dataset/ReducedLP/scpc2_reduced.lp",
"../Dataset/ReducedLP/scpc3_reduced.lp",

#"../Dataset/ReducedLP/scpclr10_reduced.lp",
#"../Dataset/ReducedLP/scpclr11_reduced.lp", #

"../Dataset/ReducedLP/scpd1_reduced.lp",
"../Dataset/ReducedLP/scpd2_reduced.lp",
"../Dataset/ReducedLP/scpd3_reduced.lp",
"../Dataset/ReducedLP/scpd4_reduced.lp",

"../Dataset/ReducedLP/scpnre1_reduced.lp",
"../Dataset/ReducedLP/scpnre2_reduced.lp",
"../Dataset/ReducedLP/scpnre3_reduced.lp",
#"../Dataset/ReducedLP/scpnre4_reduced.lp",
#"../Dataset/ReducedLP/scpnre5_reduced.lp",

"../Dataset/ReducedLP/scpnrf1_reduced.lp",
"../Dataset/ReducedLP/scpnrf2_reduced.lp",
"../Dataset/ReducedLP/scpnrf3_reduced.lp",
#"../Dataset/ReducedLP/scpnrf4_reduced.lp",
#"../Dataset/ReducedLP/scpnrf5_reduced.lp",

#"../Dataset/ReducedLP/scpnrg1_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrg2_reduced.lp",



##never used
#"../Dataset/ReducedLP/scpnrg3_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrg4_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrg5_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh2_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh3_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh5_reduced.lp", ##
#"../Dataset/ReducedLP/scpclr12_reduced.lp", ##
#"../Dataset/ReducedLP/scpclr13_reduced.lp", ##

]





#files=[
#"../Dataset/ReducedLP/scpb2_reduced.lp",
#"../Dataset/ReducedLP/scpc2_reduced.lp",
#"../Dataset/ReducedLP/scpclr10_reduced.lp",
#"../Dataset/ReducedLP/scpd2_reduced.lp",
#"../Dataset/ReducedLP/scpd4_reduced.lp",
#"../Dataset/ReducedLP/scpnre1_reduced.lp",
#"../Dataset/ReducedLP/scpnre3_reduced.lp",
#"../Dataset/ReducedLP/scpnrf3_reduced.lp",
#"../Dataset/ReducedLP/scpnrf5_reduced.lp",
#"../Dataset/ReducedLP/scpnrg2_reduced.lp",


##"../Dataset/ReducedLP/scpb4_reduced.lp",
##"../Dataset/ReducedLP/scpc3_reduced.lp",
##"../Dataset/ReducedLP/scpclr11_reduced.lp", #
##"../Dataset/ReducedLP/scpd1_reduced.lp",
##"../Dataset/ReducedLP/scpd3_reduced.lp",
##"../Dataset/ReducedLP/scpnre2_reduced.lp",
##"../Dataset/ReducedLP/scpnre4_reduced.lp",
##"../Dataset/ReducedLP/scpnre5_reduced.lp",
##"../Dataset/ReducedLP/scpnrf1_reduced.lp",
##"../Dataset/ReducedLP/scpnrf2_reduced.lp",
##"../Dataset/ReducedLP/scpnrf4_reduced.lp",
##"../Dataset/ReducedLP/scpnrg1_reduced.lp", ##

#]

#files=[
#"../Dataset/ReducedLP/scpb4_reduced.lp",
#"../Dataset/ReducedLP/scpc3_reduced.lp",
#"../Dataset/ReducedLP/scpclr11_reduced.lp", #
#"../Dataset/ReducedLP/scpd1_reduced.lp",
#"../Dataset/ReducedLP/scpd3_reduced.lp",
#"../Dataset/ReducedLP/scpnre2_reduced.lp",
#"../Dataset/ReducedLP/scpnre4_reduced.lp",
#"../Dataset/ReducedLP/scpnre5_reduced.lp",
#"../Dataset/ReducedLP/scpnrf1_reduced.lp",
#"../Dataset/ReducedLP/scpnrf2_reduced.lp",
#"../Dataset/ReducedLP/scpnrf4_reduced.lp",
#"../Dataset/ReducedLP/scpnrg1_reduced.lp", ##
#]


callback = 1

branching = 3

constraintBranchVer = 1

timelimit = 86400

#delta=[1.5]#, 2, 5]
#delta = [1.5, 2, 5]
delta = [1]
#lookAhead = [50]#, 100, 150, 200]
#lookAhead = [2147483647]#, 100, 150, 200]
lookAhead = [0]
threads = 0

#seed = [1867, 3363, 6822, 8927, 9397]

#seed = [1867, 3363, 6822, 8927, 9397, 338, 1997, 5538, 8007, 8957]

seed = [338, 1867, 6822, 8927, 9397]

#reverse = [1,]
reverse = [0]

repeated_first = [0,]

sort = [0]
#sort = [0, 1]
#sort = [0, 1, 2]



bash=open("runJobs.sh", "w")
bash.write("#!/bin/bash\n")


storeResults = 1
print(files)
counter = 0
for a in delta:
    for b in lookAhead:
        for c in seed:
            for d in reverse:
                for e in repeated_first:
                    for f in sort:
                        for g in files:
                            name="./jobs/job"+str(counter)+".slurm"
                            out = open(name, "w")
                            out.write("#!/bin/bash\n\n# specify the number of task computed together\n#SBATCH --ntasks 1\n\n# specify the number of cpus in each task\n#SBATCH --cpus-per-task 4\n\n# slurm partition try rop, arrow, razor \n#SBATCH --partition razor\n\n#SBATCH --time 30:00:00\n#SBATCH --mem 10G\n\n#SBATCH --job-name SCPconstr\n#SBATCH --output output/output_%j.txt\n#SBATCH --error errors/errors_%j.txt\n#SBATCH --mail-user peruccien@dei.unipd.it\n\ncd /home/petruccien/Constraint-branching-for-Set-Covering-Problem/build\n\n")
                            out.write("srun ./main.o -file "+ g + " -callback "+ str(callback) + " -branching " + str(branching) + " -constraintBranchVer " + str(constraintBranchVer) + " -delta " + str(a) + " -seed " + str(c) + " -time_limit "+ str(timelimit) + " -threads " + str(threads) + " -storeResults " + str(storeResults) + " -lookAhead "+ str(b) + " -reverse " + str(d) + " -repeatedFirst " + str(e) + " -sort " + str(f) + " -average " + "0\n")
                            bash.write("sbatch ../jobs/job"+str(counter)+".slurm\n")
                            counter+=1
                            if(f>0):
                                name="./jobs/job"+str(counter)+".slurm"
                                out = open(name, "w")
                                out.write("#!/bin/bash\n\n# specify the number of task computed together\n#SBATCH --ntasks 1\n\n# specify the number of cpus in each task\n#SBATCH --cpus-per-task 4\n\n# slurm partition try rop, arrow, razor \n#SBATCH --partition razor\n\n#SBATCH --time 24:00:00\n#SBATCH --mem 10G\n\n#SBATCH --job-name SCPconstr\n#SBATCH --output output/output_%j.txt\n#SBATCH --error errors/errors_%j.txt\n#SBATCH --mail-user peruccien@dei.unipd.it\n\ncd /home/petruccien/Constraint-branching-for-Set-Covering-Problem/build\n\n")
                                out.write("srun ./main.o -file "+ g + " -callback "+ str(callback) + " -branching " + str(branching) + " -constraintBranchVer " + str(constraintBranchVer) + " -delta " + str(a) + " -seed " + str(c) + " -time_limit "+ str(timelimit) + " -threads " + str(threads) + " -storeResults " + str(storeResults) + " -lookAhead "+ str(b) + " -reverse " + str(d) + " -repeatedFirst " + str(e) + " -sort " + str(f) + " -average " + "1\n")
                                bash.write("sbatch ../jobs/job"+str(counter)+".slurm\n")
                                counter+=1

#srun ./main.o -file ../Dataset/ReducedLP/scpc3_reduced.lp -callback 1 -branching 1 -constraintBranchVer 3 -delta 1.5 -seed 0 -timelimit 3600 -threads 0 -storeResults 1 -lookAhead 10
