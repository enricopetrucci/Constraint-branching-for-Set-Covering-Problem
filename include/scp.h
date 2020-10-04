/**
 * 
 * @author Petrucci Enrico
*/


#ifndef scp_H_  
#define scp_H_

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h> 
#include <stdio.h>
#include <unistd.h>
#include "chrono.h"  

#include <cplex.h>

#define EPS 1e-5	//CPLEX precision

               
/**
 * Data structure containing everything useful to an execution 
 * and the solution computed by the CPLEX library.
*/
typedef struct {   

	// execution parameters 
	double timelimit;					// overall time limit, in sec.s
	char input_file[1000];		 		// input file
	int seed;							// seed that is used by Cplex to generate random numbers in order to take decisions
	int plot;						    // if set to one it enables plots during execution
	int callback;						// 0 no callback, 1 legacy callback, 2 generic callback
	int num_cols;
	int num_rows;
	int branchcalls[2];
	int extractPreprocessing; 			// only used to extract and save the preprocessed problem at the first branching, the problem is not solved
} instance;

// functions to compute the solutions to the TSP using the CPLEX library
int scpopt(instance *inst);


#endif   /* scp_H_ */
