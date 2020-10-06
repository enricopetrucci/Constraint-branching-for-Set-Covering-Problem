/*
 * This file contains utilities functions to print errors and arrays
 * and plots using gnuplot.
 * 
 * @author Miani Eleonora
 * @author Petrucci Enrico
*/

#include "utilities.h"

/**
 * Prints a standardized error message.
 * 
 * @param err error to print
*/
void print_error(const char *err)
{
	printf("\n\n ERROR: %s \n\n", err);
	fflush(NULL);
	exit(1);
}


/**
 *  Prints the content of an array of double which has length len.
 * 
 * @param arr array to print
 * @param len length of the array to print
*/
void print_array(double *arr, int len)
{
	printf("\n[");
	for (int i=0; i < len; i++)
	{
		printf("%f, ", arr[i]);
	}
	printf("]\n");
}


/**
 *  Prints the content of an array of int which has length len.
 * 
 * @param arr array to print
 * @param len length of the array to print
*/
void print_array_int(int *arr, int len)
{
	printf("\n[");
	for (int i=0; i < len; i++)
	{
		printf("%d, ", arr[i]);
	}
	printf("]\n");
}


void fprint_array_int(FILE *f, int *arr, int len)
{
	fprintf(f,"\n[");
	for (int i=0; i < len; i++)
	{
		fprintf(f, "%d, ", arr[i]);
	}
	fprintf(f, "]\n");
}
