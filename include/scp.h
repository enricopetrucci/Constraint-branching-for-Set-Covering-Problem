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
	int callback;						// 0 generic callback, 1 legacy callback 
	int branching;						// select the kind of branching to use
	
	int threads;						// number of threads used in the computation
	int constraintBranchVer;			// version for the branching constraint
	int extractPreprocessing; 			// only used to extract and save the preprocessed problem at the first branching, the problem is not solved
	int scoreComparison; 			// only used to compare the scores of constraint branching and variable branching


	int num_cols;
	int num_rows;

	// datastructure containing the constraints intersection info
	int ** intersections;				// array that contains the pointers to arrays containing the variables in each possible intersections between constraints. It contains numIntersection pointers
	int numIntersections;				// number of non zeros inside intersections
	int * intersectionsLengths;

	double startTime;	
	int* constraintBranching;
	int* defaultBranching;

	double* timeInCallback;
	double* timeFindingConstraint;


	double* solution;					// array containing the final solution
	int shortestConstraint;
	int lowestNumVariables;
	int* variableFreq;

	int ** varConstrTable;
    int* constraintCounter; 
    int ** varConstrDim;
 

} instance;

// functions to compute the solutions to the TSP using the CPLEX library
int scpopt(instance *inst);

int preprocessinglegacycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                         const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);

void populateIntersectionsOf2(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf2Sorted(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf3(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf4(instance* inst);

void performPreprocessing(CPXENVptr env, CPXLPptr lp, instance* inst);

void findBranchingConstraint0(int* n, int* m, double* x, instance* inst);

void findBranchingConstraint1(int* n, int* m, double* x, instance* inst, double delta);

void findBranchingConstraint2(int* n, int* m, double* x, instance* inst, int* fixed);

void findBranchingConstraint3(int* n, int* m, double* x, instance* inst);

void findBranchingConstraint4(int* n, int* m, double* x, instance* inst, int* fixed);

void findBranchingConstraint5(int* n, int* m, double* x, instance* inst, int* fixed);

void computeVariableFrequency(int* izero, int* indexes, int nnz, instance* inst);

int scoreComparison(instance *inst);

int scpPreprocessing(instance *inst);

void populateIntersectionsOf2NoDup(int* izero, int* indexes, int nnz, instance* inst);

void computeConstraintsProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat,  double* productScoreConstraints, double* scoreDown, double*scoreUp, double epsilon);

void computeVariablesProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, double* productScoreVariables, double* pseudocostDown, double* pseudocostUp, double epsilon);

void computePrevisionConstraintsScores(instance* inst, double* rootSolution, double* pseudocostDown, double* pseudocostUp, double* constraintScorePrevision, double* estimateScoreDown, double* estimateScoreUp, double epsilon, int policy);

int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle);

double findBestBranchingConstraint(int *n, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);

double findBestBranchingConstraintContainingVar(int *n, int bestVar, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);


void solveUsingLegacyCallback(CPXENVptr env, CPXLPptr lp, instance* inst);

void solveUsingGenericCallback(CPXENVptr env, CPXLPptr lp, instance* inst);

void populateVariableConstraintTable(instance *inst);

void addBranchingChildsReduced(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int bestVar, int i, instance* inst);

double findBestBranchingConstraintSmart(int *n, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);

void addBranchingChilds(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, instance* inst);











#endif   /* scp_H_ */
