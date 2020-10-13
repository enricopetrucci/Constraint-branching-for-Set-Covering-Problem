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
#include <utilities.h>

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
	int callback;						// 0 no callback, 1 legacy callback, 2 generic callback
	int branching;						// select the kind of branching to use
	
	int threads;						// number of threads used in the computation
	int constraintBranchVer;			// version for the branching constraint
	int extractPreprocessing; 			// only used to extract and save the preprocessed problem at the first branching, the problem is not solved


	int num_cols;
	int num_rows;

	// datastructure containing the constraints intersection info
	int ** intersections;				// array that contains the pointers to arrays containing the variables in each possible intersections between constraints. It contains numIntersection pointers
	int * interSetLen;					// the array intersections is divided into contiguous sets that contain each intersection of the same lenghts, these length are stored in interSetLen 
	int * interSetStart;				// contains the starting position of each set.
	int numInterSet;					// number of non zeros inside interSetLen and interSetStart 
	int numIntersections;				// number of non zeros inside intersections
	
	

	double* solution;					// array containing the final solution
	int shortestConstraint;
	int lowestNumVariables;
} instance;

// functions to compute the solutions to the TSP using the CPLEX library
int scpopt(instance *inst);

int preprocessinglegacycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                         const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);

void populateIntesectionsOf2(int* izero, int* indexes, int nnz, instance* inst);

void populateIntesectionsOf3(int* izero, int* indexes, int nnz, instance* inst);

void populateIntesectionsOf4(instance* inst);

void performPreprocessing(CPXENVptr env, CPXLPptr lp, instance* inst);

void findBranchingConstraint0(int* n, int* m, double* x, instance* inst);

void findBranchingConstraint1(int* n, int* m, double* x, instance* inst, double delta);

void findBranchingConstraint2(int* n, int* m, double* x, instance* inst, int* fixed);





#endif   /* scp_H_ */
