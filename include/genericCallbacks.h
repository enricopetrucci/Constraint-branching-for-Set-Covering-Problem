/**
 * 
 * @author Petrucci Enrico
*/


#ifndef genericCallbacks_H_  
#define genericCallbacks_H_


#include "scp.h"

void findBranchingConstraint0(int* n, int* m, double* x, instance* inst);

void findBranchingConstraint1(int* n, int* m, double* x, instance* inst, double delta);

void findBranchingConstraint2(int* n, int* m, double* x, instance* inst, int* fixed);

void findBranchingConstraint3(int* n, int* m, double* x, instance* inst);

void findBranchingConstraint4(int* n, int* m, double* x, instance* inst, int* fixed);

void findBranchingConstraint5(int* n, int* m, double* x, instance* inst, int* fixed);

int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle);

void solveUsingGenericCallback(CPXENVptr env, CPXLPptr lp, instance* inst);

#endif   /* genericCallbacks_H_ */
