/**
 * This file contains the declarations of utilities functions
 * 
 * @author Petrucci Enrico
*/

#ifndef UTILITIES_H_  
#define UTILITIES_H_

#include "scp.h"

void print_error(const char *err);
void print_array(double *arr, int len);
void print_array_int(int *arr, int len);
void fprint_array_int(FILE *f, int *arr, int len);
void fprint_array_int_int(FILE *f, int ** arr, int nnz, int* izero, int* lengths, int len);



#endif   /* UTILITIES_H_ */