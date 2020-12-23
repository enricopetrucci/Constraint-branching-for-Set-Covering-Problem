/**
 * 
 * @author Petrucci Enrico
*/

#ifndef constraintIntersection_H_  
#define constraintIntersection_H_


#include "scp.h"

void populateIntersectionsOf2(int* izero, int* indexes, int nnz, instance* inst);

void computeVariableFrequencies(int *indexes, int nnz, instance *inst);

void computeConstraintScoresReducedCosts(instance *inst, double* reducedCosts);

void computeConstraintScoresFreq(instance *inst);

void populateIntersectionsOf2Original(int *izero, int *indexes, int nnz, instance *inst);

void sortIntersections(instance *inst);

void swap(instance *inst, int i, int j);

void purgeDuplicatesRepeatedFirst(instance *inst);


void purgeDuplicates(instance *inst);


int isEqual(instance *inst, int i, int j);

void merge_sort(int i, int j, int** aux, int* aux1, instance* inst);

void merge_sort1(int i, int j, int** aux, int* aux1, int* aux2,  instance* inst); 

void merge_sort2(int i, int j, int** aux, int* aux1, double * aux2,  instance* inst);

int compareIntersections(instance *inst, int i, int j);

void populateIntersectionsOf2Sorted(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf3(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf4(instance* inst);

void computeVariableFrequency(int* izero, int* indexes, int nnz, instance* inst);

void populateIntersectionsOf2NoDup(int* izero, int* indexes, int nnz, instance* inst);

void populateVariableConstraintTable(instance *inst);


#endif   /* constraintIntersection_H_ */
