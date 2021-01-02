/**
 * The purpose of this program is to solve instances of the Set Covering Problem(SCP).
 *
 * The ultimate aim is to compare variable branching with various version of constraint branching
 *   
 * @author Petrucci Enrico
*/



#include "scp.h"
#include "preprocessing.h"
#include "scoreComparison.h"



void parse_command_line(int argc, char** argv, instance* inst);



int main(int argc, char** argv)
{
	// instance of the scp problem to be solved
	instance inst;
	
	inst.startTime = second();

	// parse input from command line and populate the instance
	parse_command_line(argc, argv, &inst);

	if(inst.extractPreprocessing==1)
	{
		scpPreprocessing(&inst);   
	}
   else if(inst.scoreComparison==1)
	{
		scoreComparison(&inst);
	}
	else
	{
		// build the CPLEX model and compute the solution
		scpopt(&inst);
	}

		
	
	printf("total execution time %f sec.\n", inst.executionTime);
	return 0;
}


/**
 * Parses the arguments provided by the user when the program is called.
 * 
 * Available arguments:
 * -file, -input, -f : input file
 * -time_limit : total time limit for computation
 * -seed : random seed
 * 
 * @param argc : number of parameters
 * @param argv : pointer to the parameters
 * @param inst : instance of the TSP
*/
void parse_command_line(int argc, char** argv, instance* inst)
{
	// set the defaults: if they need to be changed it's done when parsing the arguments  
	strcpy(inst->input_file, "NULL");
	inst->timelimit = 3600;
	inst->seed = 0;
	inst->callback = 1;
	inst->extractPreprocessing = 0;
	inst->branching = 1;
	inst->constraintBranchVer = 0;
	inst->threads = 0;
	inst->scoreComparison = 0;
	inst->storeResults = 0;
	inst->delta = 2;
	inst->lookAhead = 100;
	inst->reverse = 0;
	inst->repeatedFirst = 1;
	inst->average = 1;
	inst->sort=1;


	// parse the parameters
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); printf("File: %s\n", inst->input_file); continue; } 							// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); printf("File: %s\n", inst->input_file); continue; } 							// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); printf("File: %s\n", inst->input_file); continue; } 								// input file
		if (strcmp(argv[i], "-time_limit") == 0) { inst->timelimit = atof(argv[++i]); printf("Timelimit: %f\n", inst->timelimit); continue; }						// total time limit
		if (strcmp(argv[i], "-seed") == 0) { inst->seed = atoi(argv[++i]); printf("seed: %d\n", inst->seed); continue; } 									// random seed for cplex computation
 		if (strcmp(argv[i], "-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }                      	// determine callback behavior
		if (strcmp(argv[i], "-extractPreprocessing") == 0) { inst->extractPreprocessing = atoi(argv[++i]); continue; } 	// reduce the computation to the preprocessing phase and extract the core of the problem
		if (strcmp(argv[i], "-branching") == 0) { inst->branching = atoi(argv[++i]); continue; } 						// determine branching behavior
		if (strcmp(argv[i], "-constraintBranchVer") == 0) { inst->constraintBranchVer = atoi(argv[++i]); continue; } 	// determine constraint branching behavior
		if (strcmp(argv[i], "-threads") == 0) { inst->threads = atoi(argv[++i]); continue; } 							// determine the number of threads
		if (strcmp(argv[i], "-scoreComparison") == 0) { inst->scoreComparison = atoi(argv[++i]); continue; } 			// only perform score comparison 
		if (strcmp(argv[i], "-storeResults") == 0) { inst->storeResults = atoi(argv[++i]); continue; } 					// store the results of the computation 
		if (strcmp(argv[i], "-delta") == 0) { inst->delta = atof(argv[++i]); continue; } 								// delta 
		if (strcmp(argv[i], "-lookAhead") == 0) { inst->lookAhead = atoi(argv[++i]); continue; } 					    // number of non improving constraints to be evaluated before stopping 
		if (strcmp(argv[i], "-reverse") == 0) { inst->reverse = atoi(argv[++i]); continue; } 					        // reverse = 0 constraints sorted from short to long. reverse = 1 sorted from long to short
		if (strcmp(argv[i], "-repeatedFirst") == 0) { inst->repeatedFirst = atoi(argv[++i]); continue; } 		        // if 1 the list of constraint keeps first the constraints that have duplicates 
		if (strcmp(argv[i], "-average") == 0) { inst->average = atoi(argv[++i]); continue; } 					        // average = 0 does not compute the average when for the constraint scores. 
		if (strcmp(argv[i], "-sort") == 0) { inst->sort = atoi(argv[++i]); continue; } 			  				        // sort for the score and not for the length.
		
	}
	if(inst->callback=1)
	{
		printf("Using legacy callbacks\n");
	}
	if(inst->branching!=1)
	{
		inst->constraintBranchVer = 0;
		inst->lookAhead = 0;
		inst->reverse = 0;
		inst->repeatedFirst = 0;
		inst->sort = 0;
		inst->average = 0;
		switch(inst->branching)
		{
			case 0:
				printf("CPLEX default\n");
				inst->delta = 0;
		    	break;
			case 2:
				printf("CPLEX default strong branching\n");
				inst->delta = 0;
				break;
			case 3:
				printf("Constraint branching + strong branching\n");
				inst->delta = 1;
				break;
		}
	}
	else
	{
		printf("Using branching constraint: ");
		switch(inst->constraintBranchVer)
		{
			case 0:
				printf("Using lookAhead after computing fractionlity\n");
				break;
			case 1:
				printf("all constraint considering only fractional variables\n");
				break;
			case 2:
				printf("only constraints that contains best variable\n");
				break;
			case 3:
				printf("Using lookAhead during the fractionlity computation\n");
				break;

		}
  	    if(inst->constraintBranchVer!=0 && inst->constraintBranchVer!=3)
		{
			inst->lookAhead=0;
		}
		else
		{
			printf("using lookahead = %d\n", inst->lookAhead);
		}

        if(inst->sort)
		{
			printf("sorting for score\n");	
			if(inst->average)
			{
				printf("Considering average instead of the sum for the score\n");
			}
		}
		else
		{
			printf("sorting for length\n");
			inst->average=0;
		}
		if(inst->reverse)
		{
			printf("from higher to lower\n");	
		}
		else
		{
			printf("from lower to higher\n");
		}
	}
}
