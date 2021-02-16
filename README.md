# Constraint Branching Strategies for the Set Covering Problem

This project contains the implementation of all the ideas presented in the Master's Degree thesis titled "Constraint Branching Strategies for the Set Covering Problem", written by Enrico Petrucci.
The code contained in this repository -- being this mainly a benchmarking application -- can seem confusing and repetitive at times, but in order to test little differences and changes in the implementation choices, it was required.
Other than the main benchmarking purpouse it is possible to simplify a given instance an save the result in a .LP file and evaluate Constraint Branching only at the root node by setting the arguments -extractPreprocessing and -scoreComparison to 1 respectively.
If these two arguments are not set, a benchmarking computation will start: this will lead to solving the instance given in input using a customized computation based on the CPLEX branch-and-bound.
The computation can be controlled by setting the arguments accordingly:

|Parameter|Usage|
| ------ | ------ |
|  -file <filename.LP> |  Specifies the instance to be solved | 
|  -time_limit <time limit>|  Specifies the timelimit | 
| -seed <seed> | Specifies the seed that is given to CPLEX
| -callback [0,1]| Specifies the callback type to be used. 0-> Generic Callback; 1-> Legacy Callback. All the methods considered in the thesis use Legacy Callbacks
| -branching [0, 1, 2, 3]|  0-> CPLEX default; 1-> constraint branching; 2-> CPLEX strong branching; 3-> constraint strong branching.
| -constraintBranchVer [0, 1, 2, 3, 4, 5]| Distinct version that manage differently the branching strategy
| -threads <number of threads to be used> | If 0 it automatically uses the max number of threads available
| -storeResults [0,1]| If 1 stores the all the statistics for the computation 
| -delta <delta>| Parameter that governs the choice between branching on the proposed variable instead of the best constraint
| -lookAhead <lookAhead>| Parameter that governs the look-ahead strategy behavior
| -reverse [0,1]| If set to 1 the constraints are sorted in a ascending order instead of the default descending order
| -repeatedFirst [0,1]| If set to 1 it prioritize the duplicates disjunctions
| -average [0,1] | If set to 1 uses the average score instead of simply considering the sum 
| -sort [0,1,2] | 0->length; 1-> variable frequency; 2-> reduced costs
