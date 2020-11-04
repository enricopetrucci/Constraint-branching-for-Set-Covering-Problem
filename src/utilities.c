/*
 * This file contains utilities functions
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
 *  Prints the content of an array of double which has length len.
 * 
 * @param arr array to print
 * @param len length of the array to print
*/
void print_array_char(char *arr, int len)
{
	printf("\n[");
	for (int i=0; i < len; i++)
	{
		printf("%c, ", arr[i]);
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


/**
 *  Prints the content of an array of int which has length len to file.
 * 
 * @param f file
 * @param arr array to print
 * @param len length of the array to print
*/

void fprint_array_int(FILE *f, int *arr, int len)
{
	fprintf(f,"\n[");
	for (int i=0; i < len; i++)
	{
		fprintf(f, "%d, ", arr[i]);
	}
	fprintf(f, "]\n");
}

/**
 *  Prints the content of an array of int which has length len to file.
 * 
 * @param f file
 * @param arr array to print
 * @param len length of the array to print
*/

void fprint_array(FILE *f, double *arr, int len)
{
	fprintf(f,"\n[");
	for (int i=0; i < len; i++)
	{
		fprintf(f, "%f, ", arr[i]);
	}
	fprintf(f, "]\n");
}


/**
 *  Prints the content of an array of arrays of int which has length len.
 * 
 * @param arr array to print
 * @param nnz number of elements in arr
 * @param izero array containing beginning index for each set of elements with the same lenght
 * @param lengths array containing length for each set
 * @param len length of the array to print
*/
void print_array_int_int(int ** arr, int nnz, int* izero, int* lengths, int len)
{
   for(int i=0; i < len; i++)
   {
      printf("length %d:", lengths[i]);
      for(int j = 0; j < ((i == len-1) ? nnz : izero[i+1]) - izero[i]; j++)
      {
         print_array_int(arr[izero[i]+j], lengths[i]);
      }

   }
}


/**
 *  Prints the content of an array of arrays of int which has length len.
 * 
 * @param arr array to print
 * @param nnz number of elements in arr
 * @param izero array containing beginning index for each set of elements with the same lenght
 * @param lengths array containing length for each set
 * @param len length of the array to print
*/
void print_array_int_int1(int ** arr, int** lengths, int* counter, int len)
{
   for(int i = 0; i < len; i++)
   {
	   	if(counter[i]>0)
		{
			printf("variable %d in constraints:",i);
			print_array_int(arr[i], counter[i]);
			printf("With lengths:");
			print_array_int(lengths[i], counter[i]);	
		}
   }
}


/**
 *  Prints the content of an array of arrays of int which has length len.
 * 
 * @param arr array to print
 * @param nnz number of elements in arr
 * @param izero array containing beginning index for each set of elements with the same lenght
 * @param lengths array containing length for each set
 * @param len length of the array to print
*/
void fprint_array_int_int1(FILE* f, int ** arr, int** lengths, int* counter, int len)
{
   for(int i = 0; i < len; i++)
   {
	   	if(counter[i]>0)
		{
			fprintf(f, "variable %d in constraints:",i);
			fprint_array_int(f, arr[i], counter[i]);
			fprintf(f, "With lengths:");
			fprint_array_int(f, lengths[i], counter[i]);	
		}
   }
}


/**
 *  Prints the content of an array of arrays of int which has length len to file.
 * 
 * @param f file
 * @param arr array to print
 * @param nnz number of elements in arr
 * @param izero array containing beginning index for each set of elements with the same lenght
 * @param lengths array containing length for each set
 * @param len length of the array to print
*/
void fprint_array_int_int(FILE *f, int ** arr, int nnz, int* izero, int* lengths, int len)
{
   for(int i=0; i < len; i++)
   {
      fprintf(f, "length %d:", lengths[i]);
      for(int j = 0; j < ((i == len-1) ? nnz : izero[i+1]) - izero[i]; j++)
      {
         fprint_array_int(f, arr[izero[i]+j], lengths[i]);
      }

   }
}
