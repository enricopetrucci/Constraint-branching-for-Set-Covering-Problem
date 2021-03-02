 
files=[
"../Dataset/ReducedLP/scpb2_reduced.lp",
"../Dataset/ReducedLP/scpb4_reduced.lp",

"../Dataset/ReducedLP/scpc2_reduced.lp",
"../Dataset/ReducedLP/scpc3_reduced.lp",

"../Dataset/ReducedLP/scpclr10_reduced.lp",
"../Dataset/ReducedLP/scpclr11_reduced.lp", #

"../Dataset/ReducedLP/scpd1_reduced.lp",
"../Dataset/ReducedLP/scpd2_reduced.lp",
"../Dataset/ReducedLP/scpd3_reduced.lp",
"../Dataset/ReducedLP/scpd4_reduced.lp",

"../Dataset/ReducedLP/scpnre1_reduced.lp",
"../Dataset/ReducedLP/scpnre2_reduced.lp",
"../Dataset/ReducedLP/scpnre3_reduced.lp",
"../Dataset/ReducedLP/scpnre4_reduced.lp",
"../Dataset/ReducedLP/scpnre5_reduced.lp",

"../Dataset/ReducedLP/scpnrf1_reduced.lp",
"../Dataset/ReducedLP/scpnrf2_reduced.lp",
"../Dataset/ReducedLP/scpnrf3_reduced.lp",
"../Dataset/ReducedLP/scpnrf4_reduced.lp",
"../Dataset/ReducedLP/scpnrf5_reduced.lp",

"../Dataset/ReducedLP/scpnrg1_reduced.lp", ##
"../Dataset/ReducedLP/scpnrg2_reduced.lp",


#"../Dataset/ReducedLP/scpnrg3_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrg4_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrg5_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh2_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh3_reduced.lp", ##
#"../Dataset/ReducedLP/scpnrh5_reduced.lp", ##
#"../Dataset/ReducedLP/scpclr12_reduced.lp", ##
#"../Dataset/ReducedLP/scpclr13_reduced.lp", ##

]

callback = 1

branching = 0

constraintBranchVer = 3

timelimit = 3600

threads = 0

#seed = [1867, 3363, 6822, 8927, 9397]

#seed = [338, 1997, 5538, 8007, 8957]

seed = [338, 1997, 6822, 8007, 8957]



bash=open("runJobs.sh", "w")
bash.write("#!/bin/bash\n")


storeResults = 1
print(files)
counter = 0
for c in seed:
    for g in files:
        name="./jobsDefault/job"+str(counter)+".slurm"
        out = open(name, "w")
        out.write("#!/bin/bash\n\n# specify the number of task computed together\n#SBATCH --ntasks 1\n\n# specify the number of cpus in each task\n#SBATCH --cpus-per-task 4\n\n# slurm partition try rop, arrow, razor \n#SBATCH --partition razor\n\n#SBATCH --time :00:00\n#SBATCH --mem 10G\n\n#SBATCH --job-name SCPconstr\n#SBATCH --output output/output_%j.txt\n#SBATCH --error errors/errors_%j.txt\n#SBATCH --mail-user peruccien@dei.unipd.it\n\ncd /home/petruccien/Constraint-branching-for-Set-Covering-Problem/build\n\n")
        out.write("srun ./main.o -file "+ g + " -callback "+ str(callback) + " -branching " + str(branching) + " -seed " + str(c) + " -timelimit "+ str(timelimit) + " -threads " + str(threads) + " -storeResults " + str(storeResults)+"\n")
        bash.write("sbatch ../jobsDefault/job"+str(counter)+".slurm\n")
        counter+=1

#srun ./main.o -file ../Dataset/ReducedLP/scpc3_reduced.lp -callback 1 -branching 1 -constraintBranchVer 3 -delta 1.5 -seed 0 -timelimit 3600 -threads 0 -storeResults 1 -lookAhead 10
