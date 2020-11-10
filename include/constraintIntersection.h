/**
 * 
 * @author Petrucci Enrico
*/

#ifndef constraintIntersection_H_  
#define constraintIntersection_H_


#include "scp.h"

void populateIntersectionsOf2(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf2Sorted(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf3(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf4(instance* inst);

void computeVariableFrequency(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf2NoDup(int* izero, int* indexes, int nnz, instance* inst);

void populateVariableConstraintTable(instance *inst);


#endif   /* constraintIntersection_H_ */
