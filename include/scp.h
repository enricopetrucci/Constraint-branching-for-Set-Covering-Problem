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
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


#include "chrono.h"  

#include <cplex.h>
#include <utilities.h>

#define EPS 1e-5	//CPLEX precision


/**
 * Data structure containing information about a single node during the computation.
*/

typedef struct
{
	int num;
	double value;
	int depth;
	double varFracPerc;
} nodeInfo;


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
	double endTime;
	int* constraintBranching;
	int* defaultBranching;

	double* timeInCallback;
	double* timeFindingConstraint;
	double initTime;
	double firstLoop;
	double secondLoop;

	double* solution;					// array containing the final solution
	int shortestConstraint;
	int lowestNumVariables;
	int* variableFreq;

	int ** varConstrTable;
    int* constraintCounter; 
    int ** varConstrDim;
 
	double executionTime;
	double bestInt;
	double bestVal;
	double MIPgap;
	int exploredNodes;
	int remainingNodes;
	int totalConstraintBranching;
	int totalVariableBranching;
	int storeResults;
	double delta;

	double prepTime;
	double callbackTime1;
	double callbackTime2;
	double callbackTime3;

	int* interSetLen;
	int* interSetStart;
	int numInterSet;

	// arrays to be used inside the callbacks
	double **xs;
	double **varPseudocostsUp;
	double **varPseudocostsDown;
	double **estimateScoreDown;
    double **estimateScoreUp;
    double **sum;
	double *values;
	char *lu;
	double *bd;
	int maxConstrLen;

	int* variableScores;
	int* constraintScores;

	double* constraintScoresD;

	int repeatedNum;
	int repeatedFirst;
	int reverse;
	int average;

	int sort;

	int lookAhead;

	int* constraintBranchingDepth;
	int* variableBranchingDepth;
	CPXENVptr* envs;													//CPLEX enviromnments, one per thread used for strong constraint branching
    CPXLPptr* lps;
	int** cstats;
	int** rstats;
	int* cstatDims;
	int* rstatDims;
	
	double** ogBd;
	char** ogLu;
	int** ogIndices;

	nodeInfo** nInfo;
	int* nInfoLengths;
	int* nInfoIndex;

} instance;



// functions to compute the solutions to the TSP using the CPLEX library
int scpopt(instance *inst);

double findBestBranchingConstraint(int *n, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);

double findBestBranchingConstraintContainingVar(int *n, int bestVar, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);

double findBestBranchingConstraintFracVar(int *n, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst);

double findBestBranchingConstraintFracVarLookAhead(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst, int threadNum);

double findBestBranchingConstraintLookAhead(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst);

double findBestBranchingConstraintStrongContVar(CPXCENVptr env, void *cbdata, int wherefrom, int *n, int bestVar, int threadNum, double obj, instance *inst);

void solveUsingLegacyCallback(CPXENVptr env, CPXLPptr lp, instance* inst);

void addBranchingChilds(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, int threadNum, instance *inst);

void addBranchingChildsAlloc(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, instance *inst);


void saveComputationResults(instance *inst);


int plotTreeCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                            const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);


void resetPlotTree(instance *inst);

int legacyBranchingCallbackStrong(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                            const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p);

void prepareBranchingConstraints(CPXCENVptr env, CPXLPptr lp, instance *inst);

void sortIntersectionWRTreducedCosts(CPXCENVptr env, void *cbdata, int wherefrom, instance *inst);

double computeStrongBranchingOnConstraint(instance* inst, double obj, int* cstat, int* rstat, int ogNumCols, int threadNum, int i);

void getBase(CPXCENVptr env, CPXLPptr lp, instance* inst, int threadNum, int cur_numcols, int cur_numrows);

double computeVariableProductScore(instance* inst, double* rootSolution, double obj, int threadNum, int bdcnt, const int *nodebeg, const int* indices, const char *lu, const double *bd);

void saveNodeInfoToFile(instance *inst);


#endif   /* scp_H_ */
