/**
 * The purpose of this program is to solve instances of the Set Covering Problem(SCP).
 *
 * The ultimate aim is to compare variable branching with various version of constraint branching
 *   
 * @author Petrucci Enrico
*/



#include "scp.h"

void parse_command_line(int argc, char** argv, instance* inst);



int main(int argc, char** argv)
{
	// instance of the scp problem to be solved
	instance inst;

	// parse input from command line and populate the instance
	parse_command_line(argc, argv, &inst);


		
	// build the CPLEX model and compute the solution
	scpopt(&inst);
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
	inst->timelimit = 100000;
	inst->seed = 0;
	inst->plot = 0;
	inst->callback = 1;
	inst->extractPreprocessing = 0;

	// parse the parameters
	for (int i = 1; i < argc; i++)
	{
		if (strcmp(argv[i], "-file") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 					// input file
		if (strcmp(argv[i], "-input") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 					// input file
		if (strcmp(argv[i], "-f") == 0) { strcpy(inst->input_file, argv[++i]); continue; } 						// input file
		if (strcmp(argv[i], "-time_limit") == 0) { inst->timelimit = atof(argv[++i]); continue; }				// total time limit
		if (strcmp(argv[i], "-seed") == 0) { inst->seed = atoi(argv[++i]); continue; } 							// random seed for cplex computation
 		if (strcmp(argv[i], "-plot") == 0) { inst->plot = atoi(argv[++i]); continue; }                  		// plot during execution
		if (strcmp(argv[i], "-callback") == 0) { inst->callback = atoi(argv[++i]); continue; }                      // determine callback behavior
		if (strcmp(argv[i], "-extractPreprocessing") == 0) { inst->extractPreprocessing = atoi(argv[++i]); continue; } // reduce the computation to the preprocessing phase and extract the core of the problem

	}
}
