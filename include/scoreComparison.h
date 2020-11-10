/**
 * 
 * @author Petrucci Enrico
*/


#ifndef scoreComparison_H_  
#define scoreComparison_H_

#include "scp.h"

int scoreComparison(instance *inst);

void computeConstraintsProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat,  double* productScoreConstraints, double* scoreDown, double*scoreUp, double epsilon);

void computeVariablesProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, double* productScoreVariables, double* pseudocostDown, double* pseudocostUp, double epsilon);

void computePrevisionConstraintsScores(instance* inst, double* rootSolution, double* pseudocostDown, double* pseudocostUp, double* constraintScorePrevision, double* estimateScoreDown, double* estimateScoreUp, double epsilon, int policy);

#endif   /* scoreComparison_H_ */
