/**
 * This file contains functions needed for managing the constraint intersections and other related datastructures
 * scp using the CPLEX library.
 *
 * @author Petrucci Enrico
*/

#include "scp.h"
#include "constraintIntersection.h"

/**
 * Compute all possible intersections between the
 * constraints taken in paris.
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void populateIntersectionsOf2(int *izero, int *indexes, int nnz, instance *inst)
{
    int **intersections = (int **)calloc(inst->num_rows * (inst->num_rows - 1) / 2, sizeof(int *));
    int *intersectionsLengths = (int *)calloc(inst->num_rows * (inst->num_rows - 1) / 2, sizeof(int));
    // initialize estimateScoreUp to the max int value

    int numIntersections = 0;

    // cycle the constraints, choosing the first one
    for (int i = 0; i < inst->num_rows - 1; i++)
    {
        //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
        //Select the second constraint
        for (int n = i + 1; n < inst->num_rows; n++)
        {
            int j = 0;
            int m = 0;
            int intersectionLen = 0;

            int length = (n < inst->num_rows - 1 ? izero[n + 1] - izero[n] : nnz - izero[n]);
            //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

            // alloc space for current intersection
            intersections[numIntersections] = (int *)calloc(inst->num_cols, sizeof(int));

            // cycle on all the variables inside both constraints
            while (j < izero[i + 1] - izero[i] && m < length) // problem for the last constraint's length
            {
                // case 1: same variable. The variable is added to the computed information
                if (indexes[izero[i] + j] == indexes[izero[n] + m])
                {
                    // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                    intersections[numIntersections][intersectionLen] = indexes[izero[i] + j];
                    intersectionLen++;
                    m++;
                    j++;
                }
                else
                {
                    // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
                    if (indexes[izero[i] + j] < indexes[izero[n] + m])
                    {
                        j++;
                    }
                    else
                    {
                        if (indexes[izero[n] + m] < indexes[izero[i] + j])
                        {
                            m++;
                        }
                    }
                }
            }
            // if the intersecion is non empty resize array
            if (intersectionLen > 1)
            {
                
                intersections[numIntersections] = (int *)realloc(intersections[numIntersections], intersectionLen * sizeof(int));
                intersectionsLengths[numIntersections] = intersectionLen;
                numIntersections++;
            }
            else
            {
                free(intersections[numIntersections]);
            }
        }
    }
    //print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);

    //printf("Total length %d max length %d\n", numIntersections, inst->num_rows * (inst->num_rows - 1) / 2);

    inst->intersections = intersections;
    inst->intersectionsLengths = intersectionsLengths;
    inst->numIntersections = numIntersections;
    
    // FILE *f;
    // f = fopen("IntersectionsUnsorted.txt", "w");
    // fprint_array_int_int2(f, intersections, intersectionsLengths, numIntersections);
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
}

// function to sort the subsection a[i .. j] of the array a[]
void merge_sort(int i, int j, int** aux, instance* inst) {
    int** a = inst->intersections;
    if (j <= i) {
        return;     // the subsection is empty or a single element
    }
    int mid = (i + j) / 2;

    // left sub-array is a[i .. mid]
    // right sub-array is a[mid + 1 .. j]
    merge_sort(i, mid, aux, inst);     // sort the left sub-array recursively
    merge_sort(mid + 1, j, aux, inst);     // sort the right sub-array recursively
    
    int pointer_left = i;       // pointer_left points to the beginning of the left sub-array
    int pointer_right = mid + 1;        // pointer_right points to the beginning of the right sub-array
    int k;      // k is the loop counter

    // we loop from i to j to fill each element of the final merged array
    for (k = i; k <= j; k++) {
        if (pointer_left == mid + 1) {      // left pointer has reached the limit
            aux[k] = a[pointer_right];
            pointer_right++;
        } else if (pointer_right == j + 1) {        // right pointer has reached the limit
            aux[k] = a[pointer_left];
            pointer_left++;
        } else if (compareIntersections(inst, pointer_right, pointer_left)) {        // pointer left points to smaller element
            aux[k] = a[pointer_left];
            pointer_left++;
        } else {        // pointer right points to smaller element
            aux[k] = a[pointer_right];
            pointer_right++;
        }
    }

    for (k = i; k <= j; k++) {      // copy the elements from aux[] to a[]
        a[k] = aux[k];
    }
    // FILE *f;
    // f = fopen("sortedConstraints.txt", "w");
    // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
    // printf("breakpoint\n");
}



// function to sort the subsection a[i .. j] of the array a[]
void merge_sort1(int i, int j, int** aux, int* aux1, instance* inst) {
    int** a = inst->intersections;
    int* a1 = inst->intersectionsLengths;

    if (j <= i) {
        return;     // the subsection is empty or a single element
    }
    int mid = (i + j) / 2;

    // left sub-array is a[i .. mid]
    // right sub-array is a[mid + 1 .. j]
    merge_sort1(i, mid, aux, aux1,  inst);     // sort the left sub-array recursively
    merge_sort1(mid + 1, j, aux, aux1, inst);     // sort the right sub-array recursively
    
    int pointer_left = i;       // pointer_left points to the beginning of the left sub-array
    int pointer_right = mid + 1;        // pointer_right points to the beginning of the right sub-array
    int k;      // k is the loop counter

    // we loop from i to j to fill each element of the final merged array
    for (k = i; k <= j; k++) {
        if (pointer_left == mid + 1) {      // left pointer has reached the limit
            aux[k] = a[pointer_right];
            aux1[k]=a1[pointer_right];
            pointer_right++;
        } else if (pointer_right == j + 1) {        // right pointer has reached the limit
            aux[k] = a[pointer_left];
            aux1[k] = a1[pointer_left];
            pointer_left++;
        } else if (compareIntersections(inst, pointer_right, pointer_left)) {        // pointer left points to smaller element
            aux[k] = a[pointer_left];
            aux1[k] = a1[pointer_left];
            pointer_left++;
        } else {        // pointer right points to smaller element
            aux[k] = a[pointer_right];
            aux1[k] = a1[pointer_right];
            pointer_right++;
        }
    }

    for (k = i; k <= j; k++) {      // copy the elements from aux[] to a[]
        a[k] = aux[k];
        a1[k] = aux1[k];
    }
    // FILE *f;
    // f = fopen("sortedConstraints.txt", "w");
    // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
    // printf("breakpoint\n");
}


void sortIntersections(instance *inst)
{
    int i, j; 
    int min_idx;

    for (i = 0; i < inst->numIntersections - 1; i++)
    { 
        min_idx = i;
        for(j = i; j < inst->numIntersections; j++)
        {
            if(compareIntersections(inst, min_idx, j))
            {
                min_idx = j;
            }
        }
        swap(inst, i, min_idx);
    }
}


void swap(instance *inst, int i, int j)
{
    int* tempInter = inst->intersections[i];
    int tempLen = inst->intersectionsLengths[i];
    
    inst->intersections[i]=inst->intersections[j];
    inst->intersectionsLengths[i]=inst->intersectionsLengths[j];

    inst->intersections[j]=tempInter;
    inst->intersectionsLengths[j]=tempLen;
}


int compareIntersections(instance *inst, int i, int j)
{
    // j<i for length
    if(inst->intersectionsLengths[i] > inst->intersectionsLengths[j])
        return 1;
    if(inst->intersectionsLengths[i] < inst->intersectionsLengths[j])
        return 0;
    for (int k=0; k<inst->intersectionsLengths[i]; k++)
    {
        if(inst->intersections[i][k] > inst->intersections[j][k])
            return 1;
        if(inst->intersections[i][k] < inst->intersections[j][k])
            return 0;
    }
    return 0;
}



void purgeDuplicates(instance *inst)
{
    int **intersections = (int **)calloc(inst->numIntersections, sizeof(int *));
    int *intersectionsLen =(int*)calloc(inst->numIntersections, sizeof(int));
    
    int lastIn = 0; // index in the old array of the last constraint inserted.
    int idxNew = 0; // index in the new array in which we want to add the next entry 
    int nDup = 0;
    
    int repeatedIdx = 0;
    int repeating = 1;

    intersections[idxNew] = inst->intersections[0]; 
    intersectionsLen[idxNew++] = inst->intersectionsLengths[0];
    for(int i=1; i<inst->numIntersections; i++)
    {
        while(isEqual(inst, lastIn, i)&& i<inst->numIntersections)
        {
            if(repeating)
            {
                int* temp = intersections[repeatedIdx];
                int tempLen = intersectionsLen[repeatedIdx];
                intersections[repeatedIdx]=intersections[idxNew-1];
                intersections[idxNew-1]=temp;
                intersectionsLen[repeatedIdx]=intersectionsLen[idxNew-1];
                intersectionsLen[idxNew-1]=tempLen;
                repeatedIdx++;
                repeating=0;
            }            
            free(inst->intersections[i]);
            nDup++;
            i++;
        }
        repeating=1;
        intersections[idxNew] = inst->intersections[i]; 
        intersectionsLen[idxNew++] = inst->intersectionsLengths[i];
        lastIn=i;
    
    }
    printf("Found %d duplicates\n",nDup);
    free(inst->intersections);
    free(inst->intersectionsLengths);

    inst->intersectionsLengths = intersectionsLen;
    inst->intersections = intersections;
    inst->numIntersections = idxNew;

}

int isEqual(instance *inst, int i, int j)
{
    if(inst->intersectionsLengths[i] != inst->intersectionsLengths[j])
        return 0;
    
    for (int k=0; k<inst->intersectionsLengths[i]; k++)
    {
        if(inst->intersections[i][k] != inst->intersections[j][k])
            return 0;
    }
    return 1;
}

/**
 * Populate a matrix that at each variable associate all the possible intersections of constraints in which it appears.
 * The result is stored inside inst->varConstraintTable
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void populateVariableConstraintTable(instance *inst)
{
    // 
    int **varConstrTable = (int **)calloc(inst->num_cols, sizeof(int *));
    // tells how long each constraint is 
    int *constraintCounter = (int *)calloc(inst->num_cols, sizeof(int));

    for (int i = 0; i < inst->num_cols; i++)
    {
        varConstrTable[i] = (int *)calloc(inst->numIntersections, sizeof(int));
    }

    // cycle on the constraint
    for (int i = 0; i < inst->numIntersections; i++)
    {
        // cycle on each variable inside the constraint
        for (int k = 0; k < inst->intersectionsLengths[i]; k++)
        {
            int currentVar = (inst->intersections[i][k]);
            varConstrTable[currentVar][constraintCounter[currentVar]] = i;
            constraintCounter[currentVar]++;
        }
    }
    // FILE *f;
    // f = fopen("varToConstraints.txt", "w");
    // fprint_array_int_int2(f, varConstrTable, constraintCounter, inst->num_cols);
    // fclose(f);
    inst->varConstrTable = varConstrTable;
    inst->constraintCounter = constraintCounter;
}




/**
 * Compute all possible intersections between the
 * constraint and stored them sorted by length inside inst->intersections
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void populateIntersectionsOf2Sorted(int *izero, int *indexes, int nnz, instance *inst)
{
    inst->intersections = (int **)calloc(inst->num_rows * (inst->num_rows - 1) / 2, sizeof(int *));
    int* interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
    int* interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
    int numInterSet = 0;
    inst->numIntersections = 0;

    // cycle the constraints, choosing the first one
    for (int i = 0; i < inst->num_rows - 1; i++)
    {
        //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
        //Select the second constraint
        for (int n = i + 1; n < inst->num_rows; n++)
        {
            int j = 0;
            int m = 0;
            int intersectionLen = 0;

            int length = (n < inst->num_rows - 1 ? izero[n + 1] - izero[n] : nnz - izero[n]);
            //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

            // alloc space for current intersection
            inst->intersections[inst->numIntersections] = (int *)calloc(inst->num_cols, sizeof(int));

            // cycle on all the variables inside both constraints
            while (j < izero[i + 1] - izero[i] && m < length) // problem for the last constraint's length
            {
                // case 1: same variable. The variable is added to the computed information
                if (indexes[izero[i] + j] == indexes[izero[n] + m])
                {
                    // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                    inst->intersections[inst->numIntersections][intersectionLen] = indexes[izero[i] + j];
                    intersectionLen++;
                    m++;
                    j++;
                }
                // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
                if (indexes[izero[i] + j] < indexes[izero[n] + m])
                {
                    j++;
                }

                if (indexes[izero[n] + m] < indexes[izero[i] + j])
                {
                    m++;
                }
            }
            // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
            if (intersectionLen > 1)
            {
                // realloc should free the unused spaced of the array
                inst->intersections[inst->numIntersections] = (int *)realloc(inst->intersections[inst->numIntersections], intersectionLen * sizeof(int));
                // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
                // case 1 validIntersection is empty
                if (inst->numIntersections == 0)
                {
                    // printf("Case 1 empty list\n");
                    interSetLen[0] = intersectionLen;
                    interSetStart[0] = 0;
                    numInterSet++;
                }
                else
                {
                    // cycle on all sets
                    for (int k = 0; k < numInterSet; k++)
                    {
                        // Case 2 Add new set at the beginning or add new set in between two sets
                        if ((k == 0 && intersectionLen > interSetLen[k]) || (k > 0 && intersectionLen < interSetLen[k - 1] && intersectionLen > interSetLen[k]))
                        {
                            // printf("Case 2 new set in between\n");

                            // Update information on interSetLen and interSetStart
                            for (int a = numInterSet; a >= k; a--)
                            {
                                interSetLen[a + 1] = interSetLen[a];
                                interSetStart[a + 1] = interSetStart[a] + 1;
                            }

                            interSetLen[k] = intersectionLen;
                            numInterSet++;

                            // Updae infromation in validIntersection on the right with respect to k
                            for (int a = k; a < numInterSet - 1; a++)
                            {
                                // printf("from last position to %d\n", interSetStart[a+1]-1);
                                int *temp = inst->intersections[interSetStart[a + 1] - 1];
                                inst->intersections[interSetStart[a + 1] - 1] = inst->intersections[inst->numIntersections];
                                inst->intersections[inst->numIntersections] = temp;
                            }
                            break;
                        }

                        // Case 3 no new set needed, only add an intersection and update everithing accordingly
                        if (intersectionLen == interSetLen[k])
                        {

                            // printf("Case 3 same set\n");
                            for (int a = numInterSet; a > k; a--)
                            {
                                interSetStart[a] = interSetStart[a] + 1;
                            }

                            for (int a = k + 1; a < numInterSet; a++)
                            {
                                int *temp = inst->intersections[interSetStart[a] - 1];
                                // printf("interSetStart[a]-1 = %d",interSetStart[a]-1);
                                inst->intersections[interSetStart[a] - 1] = inst->intersections[inst->numIntersections];
                                inst->intersections[inst->numIntersections] = temp;
                            }
                            break;
                        }

                        // Case 4 Add new set at the end
                        if (k == numInterSet - 1 && intersectionLen < interSetLen[k])
                        {

                            //printf("Case 4 new set at the end\n");
                            interSetLen[numInterSet] = intersectionLen;
                            interSetStart[numInterSet] = inst->numIntersections;
                            numInterSet++;
                            break;
                        }
                    }
                }
                inst->numIntersections++;
                // printf("Updated len and start:\n");
                // print_array_int(interSetLen, numInterSet);
                // print_array_int(interSetStart, numInterSet);
                // print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);
                // printf("Total length %d\n", inst->numIntersections);
                // printf("breakpoint\n");
            }
            else
            { // if the intersection is empty: free allocated space
                //printf("free\n");
                free(inst->intersections[inst->numIntersections]);
            }
        }
    }
    inst->intersections = (int **)realloc(inst->intersections, inst->numIntersections * sizeof(int *));
    inst->shortestConstraint = interSetLen[0];
    inst->lowestNumVariables = interSetLen[0];

    int *intersectionsLengths = (int *)calloc(inst->numIntersections, sizeof(int));
    // initialize estimateScoreUp to the max int value

    int set = 0;
    for (int i = 0; i < inst->numIntersections; i++)
    {
        if (set < numInterSet - 1 ? i < interSetStart[set + 1] : i < inst->numIntersections)
        {
            ;
        }
        else
            set++;

        intersectionsLengths[i] = interSetLen[set];
    }
    inst->intersectionsLengths = intersectionsLengths;
    inst->interSetLen=interSetLen;
    inst->interSetStart=interSetStart;
    inst->numInterSet=numInterSet;

    // printf("Length for each set");
    // print_array_int(interSetLen, numInterSet);
    // printf("Starting index for each set for each set");
    // print_array_int(interSetStart, numInterSet);
    // //print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);

    // printf("Total length %d max length %d\n", inst->numIntersections, inst->num_rows * (inst->num_rows - 1) / 2);

    // FILE *f;
    // f = fopen("Intersections.txt", "w");

    // fprint_array_int_int(f, inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);
    
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
}



/**
 * Compute all possible intersections between the
 * constraint and stored them sorted by length inside inst->intersections
 * without duplicates
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void populateIntersectionsOf2NoDup(int *izero, int *indexes, int nnz, instance *inst)
{
    inst->intersections = (int **)calloc(inst->num_rows * (inst->num_rows - 1) / 2, sizeof(int *));
    int* interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
    int* interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
    int numInterSet = 0;
    inst->numIntersections = 0;

    // cycle the constraints, choosing the first one
    for (int i = 0; i < inst->num_rows - 1; i++)
    {
        //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
        //Select the second constraint
        for (int n = i + 1; n < inst->num_rows; n++)
        {
            int j = 0;
            int m = 0;
            int intersectionLen = 0;

            int length = (n < inst->num_rows - 1 ? izero[n + 1] - izero[n] : nnz - izero[n]);
            //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

            // alloc space for current intersection
            inst->intersections[inst->numIntersections] = (int *)calloc(inst->num_cols, sizeof(int));

            // cycle on all the variables inside both constraints
            while (j < izero[i + 1] - izero[i] && m < length) // problem for the last constraint's length
            {
                // case 1: same variable. The variable is added to the computed information
                if (indexes[izero[i] + j] == indexes[izero[n] + m])
                {
                    // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                    inst->intersections[inst->numIntersections][intersectionLen] = indexes[izero[i] + j];
                    intersectionLen++;
                    m++;
                    j++;
                }
                // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
                if (indexes[izero[i] + j] < indexes[izero[n] + m])
                {
                    j++;
                }

                if (indexes[izero[n] + m] < indexes[izero[i] + j])
                {
                    m++;
                }
            }
            // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
            if (intersectionLen > 1)
            {
                // realloc should free the unused spaced of the array
                inst->intersections[inst->numIntersections] = (int *)realloc(inst->intersections[inst->numIntersections], intersectionLen * sizeof(int));
                // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
                // case 1 validIntersection is empty
                if (inst->numIntersections == 0)
                {
                    // printf("Case 1 empty list\n");
                    interSetLen[0] = intersectionLen;
                    interSetStart[0] = 0;
                    numInterSet++;
                    inst->numIntersections++;
                }
                else
                {
                    //check if it is already present
                    int alreadyIn = 0;
                    for (int k = 0; k < numInterSet; k++)
                    {
                        // if a subset has the same length
                        if (interSetLen[k] == intersectionLen)
                        {
                            // for each intersection in the subset
                            for (int a = interSetStart[k]; a < (k == numInterSet - 1 ? inst->numIntersections : interSetStart[k + 1]); a++)
                            {
                                int counter = 0;
                                //for each variable in a single subset
                                for (int b = 0; b < interSetLen[k]; b++)
                                {
                                    if (inst->intersections[a][b] != inst->intersections[inst->numIntersections][b])
                                    {
                                        break;
                                    }
                                    else
                                    {
                                        counter++;
                                    }
                                }
                                if (counter == interSetLen[k])
                                {
                                    alreadyIn = 1;
                                    break;
                                }
                            }
                            break;
                        }
                    }
                    if (alreadyIn == 0)
                    {
                        // cycle on all sets
                        for (int k = 0; k < numInterSet; k++)
                        {
                            // Case 2 Add new set at the beginning or add new set in between two sets
                            if ((k == 0 && intersectionLen > interSetLen[k]) || (k > 0 && intersectionLen < interSetLen[k - 1] && intersectionLen > interSetLen[k]))
                            {
                                // printf("Case 2 new set in between\n");

                                // Update information on interSetLen and interSetStart
                                for (int a = numInterSet; a >= k; a--)
                                {
                                    interSetLen[a + 1] = interSetLen[a];
                                    interSetStart[a + 1] = interSetStart[a] + 1;
                                }

                                interSetLen[k] = intersectionLen;
                                numInterSet++;

                                // Updae infromation in validIntersection on the right with respect to k
                                for (int a = k; a < numInterSet - 1; a++)
                                {
                                    // printf("from last position to %d\n", interSetStart[a+1]-1);
                                    int *temp = inst->intersections[interSetStart[a + 1] - 1];
                                    inst->intersections[interSetStart[a + 1] - 1] = inst->intersections[inst->numIntersections];
                                    inst->intersections[inst->numIntersections] = temp;
                                }
                                break;
                            }

                            // Case 3 no new set needed, only add an intersection and update everithing accordingly
                            if (intersectionLen == interSetLen[k])
                            {

                                // printf("Case 3 same set\n");
                                for (int a = numInterSet; a > k; a--)
                                {
                                    interSetStart[a] = interSetStart[a] + 1;
                                }

                                for (int a = k + 1; a < numInterSet; a++)
                                {
                                    int *temp = inst->intersections[interSetStart[a] - 1];
                                    // printf("interSetStart[a]-1 = %d",interSetStart[a]-1);
                                    inst->intersections[interSetStart[a] - 1] = inst->intersections[inst->numIntersections];
                                    inst->intersections[inst->numIntersections] = temp;
                                }
                                break;
                            }

                            // Case 4 Add new set at the end
                            if (k == numInterSet - 1 && intersectionLen < interSetLen[k])
                            {

                                //printf("Case 4 new set at the end\n");
                                interSetLen[numInterSet] = intersectionLen;
                                interSetStart[numInterSet] = inst->numIntersections;
                                numInterSet++;
                                break;
                            }
                        }
                        inst->numIntersections++;
                    }
                }

                // printf("Updated len and start:\n");
                // print_array_int(interSetLen, numInterSet);
                // print_array_int(interSetStart, numInterSet);
                // print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);
                // printf("Total length %d\n", inst->numIntersections);
                // printf("breakpoint\n");
            }
            else
            { // if the intersection is empty: free allocated space
                //printf("free\n");
                free(inst->intersections[inst->numIntersections]);
            }
        }
    }
    
    inst->intersections = (int **)realloc(inst->intersections, inst->numIntersections * sizeof(int *));
    inst->shortestConstraint = interSetLen[0];
    inst->lowestNumVariables = interSetLen[0];

    int *intersectionsLengths = (int *)calloc(inst->numIntersections, sizeof(int));
    // initialize estimateScoreUp to the max int value

    int set = 0;
    for (int i = 0; i < inst->numIntersections; i++)
    {
        if (set < numInterSet - 1 ? i < interSetStart[set + 1] : i < inst->numIntersections)
        {
            ;
        }
        else
            set++;

        intersectionsLengths[i] = interSetLen[set];
    }
    inst->intersectionsLengths = intersectionsLengths;
    // printf("Length for each set");
    // print_array_int(interSetLen, numInterSet);
    // printf("Starting index for each set for each set");
    // print_array_int(interSetStart, numInterSet);
    // //print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);

    // printf("Total length %d max length %d\n", inst->numIntersections, inst->num_rows * (inst->num_rows - 1) / 2);

    // FILE *f;
    // f = fopen("Intersections.txt", "w");

    // fprint_array_int_int(f, inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);
    
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
}

/**
 * Compute all possible intersections between the
 * constraint and stored them sorted by length inside inst->intersections
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void populateIntersectionsOf3(int *izero, int *indexes, int nnz, instance *inst)
{
    int possibleIntersections = (inst->num_rows) * (inst->num_rows - 1) * (inst->num_rows - 2) / 6;
    inst->intersections = (int **)calloc(possibleIntersections, sizeof(int *));
    int* interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
    int* interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
    int numInterSet = 0;
    inst->numIntersections = 0;

    // cycle the constraints, choosing the first one
    for (int x = 0; x < inst->num_rows - 2; x++)
    {
        //printf("x = %d\n", x);
        for (int i = x + 1; i < inst->num_rows - 1; i++)
        {
            //printf("i = %d\n",i);
            //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
            //Select the second constraint
            for (int n = i + 1; n < inst->num_rows; n++)
            {
                //printf("n = %d\n",n);
                int y = 0;
                int j = 0;
                int m = 0;
                int intersectionLen = 0;

                int length = (n < inst->num_rows - 1 ? izero[n + 1] - izero[n] : nnz - izero[n]);
                //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

                // alloc space for current intersection
                inst->intersections[inst->numIntersections] = (int *)calloc(inst->num_cols, sizeof(int));

                // cycle on all the variables inside both constraints
                while (y < izero[x + 1] - izero[x] && j < izero[i + 1] - izero[i] && m < length) // problem for the last constraint's length
                {
                    //printf("y = %d, j = %d, m = %d\n", y, j, m);
                    //printf("first variable = %d, second variable = %d, third variable = %d\n", indexes[izero[x]+y], indexes[izero[i]+j], indexes[izero[n]+m]);
                    // case 1: same variable. The variable is added to the computed information
                    if (indexes[izero[x] + y] == indexes[izero[i] + j] && indexes[izero[i] + j] == indexes[izero[n] + m])
                    {
                        //printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                        inst->intersections[inst->numIntersections][intersectionLen] = indexes[izero[i] + j];
                        intersectionLen++;
                        y++;
                        m++;
                        j++;
                    }
                    if (indexes[izero[x] + y] < indexes[izero[i] + j] && indexes[izero[x] + y] < indexes[izero[n] + m])
                    {
                        y++;
                    }
                    if (indexes[izero[i] + j] < indexes[izero[x] + y] && indexes[izero[i] + j] < indexes[izero[n] + m])
                    {
                        j++;
                    }
                    if (indexes[izero[n] + m] < indexes[izero[x] + y] && indexes[izero[n] + m] < indexes[izero[i] + j])
                    {
                        m++;
                    }
                    if (indexes[izero[x] + y] > indexes[izero[i] + j] && indexes[izero[x] + y] > indexes[izero[n] + m])
                    {
                        j++;
                        m++;
                    }
                    if (indexes[izero[i] + j] > indexes[izero[x] + y] && indexes[izero[i] + j] > indexes[izero[n] + m])
                    {
                        y++;
                        m++;
                    }
                    if (indexes[izero[n] + m] > indexes[izero[x] + y] && indexes[izero[n] + m] > indexes[izero[i] + j])
                    {
                        y++;
                        j++;
                    }
                }
                // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
                if (intersectionLen > 1)
                {
                    // realloc should free the unused spaced of the array
                    inst->intersections[inst->numIntersections] = (int *)realloc(inst->intersections[inst->numIntersections], intersectionLen * sizeof(int));
                    //printf("Intersection between constraints %d, %d and %d has %d variables\n", x, i, n, intersectionLen);
                    // case 1 validIntersection is empty
                    if (inst->numIntersections == 0)
                    {
                        // printf("Case 1 empty list\n");
                        interSetLen[0] = intersectionLen;
                        interSetStart[0] = 0;
                        numInterSet++;
                    }
                    else
                    {
                        // cycle on all sets
                        for (int k = 0; k < numInterSet; k++)
                        {
                            // Case 2 Add new set at the beginning or add new set in between two sets
                            if ((k == 0 && intersectionLen > interSetLen[k]) || (k > 0 && intersectionLen < interSetLen[k - 1] && intersectionLen > interSetLen[k]))
                            {
                                // printf("Case 2 new set in between\n");

                                // Update information on interSetLen and interSetStart
                                for (int a = numInterSet; a >= k; a--)
                                {
                                    interSetLen[a + 1] = interSetLen[a];
                                    interSetStart[a + 1] = interSetStart[a] + 1;
                                }

                                interSetLen[k] = intersectionLen;
                                numInterSet++;

                                // Update infromation in validIntersection on the right with respect to k
                                for (int a = k; a < numInterSet - 1; a++)
                                {
                                    // printf("from last position to %d\n", interSetStart[a+1]-1);
                                    int *temp = inst->intersections[interSetStart[a + 1] - 1];
                                    inst->intersections[interSetStart[a + 1] - 1] = inst->intersections[inst->numIntersections];
                                    inst->intersections[inst->numIntersections] = temp;
                                }
                                break;
                            }

                            // Case 3 no new set needed, only add an intersection and update everithing accordingly
                            if (intersectionLen == interSetLen[k])
                            {

                                // printf("Case 3 same set\n");
                                for (int a = numInterSet; a > k; a--)
                                {
                                    interSetStart[a] = interSetStart[a] + 1;
                                }

                                for (int a = k + 1; a < numInterSet; a++)
                                {
                                    int *temp = inst->intersections[interSetStart[a] - 1];
                                    // printf("interSetStart[a]-1 = %d",interSetStart[a]-1);
                                    inst->intersections[interSetStart[a] - 1] = inst->intersections[inst->numIntersections];
                                    inst->intersections[inst->numIntersections] = temp;
                                }
                                break;
                            }

                            // Case 4 Add new set at the end
                            if (k == numInterSet - 1 && intersectionLen < interSetLen[k])
                            {

                                //printf("Case 4 new set at the end\n");
                                interSetLen[numInterSet] = intersectionLen;
                                interSetStart[numInterSet] = inst->numIntersections;
                                numInterSet++;
                                break;
                            }
                        }
                    }
                    inst->numIntersections++;
                    // printf("Updated len and start:\n");
                    // print_array_int(interSetLen, numInterSet);
                    // print_array_int(interSetStart, numInterSet);
                    // print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, interSetLen, numInterSet);
                    // printf("Total length %d\n", inst->numIntersections);
                    // printf("breakpoint\n");
                }
                else
                { // if the intersection is empty: free allocated space
                    //printf("free\n");
                    free(inst->intersections[inst->numIntersections]);
                }
            }
        }
    }
    inst->intersections = (int **)realloc(inst->intersections, inst->numIntersections * sizeof(int *));
    inst->shortestConstraint = interSetLen[0];
    inst->lowestNumVariables = interSetLen[0];

    int *intersectionsLengths = (int *)calloc(inst->numIntersections, sizeof(int));
    // initialize estimateScoreUp to the max int value

    int set = 0;
    for (int i = 0; i < inst->numIntersections; i++)
    {
        if (set < numInterSet - 1 ? i < interSetStart[set + 1] : i < inst->numIntersections)
        {
            ;
        }
        else
            set++;

        intersectionsLengths[i] = interSetLen[set];
    }
    inst->intersectionsLengths = intersectionsLengths;
    // printf("Length for each set");
    // print_array_int(interSetLen, numInterSet);
    // printf("Starting index for each set for each set");
    // print_array_int(interSetStart, numInterSet);
    // //print_array_int_int(inst->intersections, inst->numIntersections, interSetStart, inst->interSetLen, numInterSet);

    // printf("Total length %d max length %d\n", inst->numIntersections, inst->num_rows * (inst->num_rows - 1) / 2);

    // FILE *f;
    // f = fopen("Intersections.txt", "w");

    // fprint_array_int_int(f, inst->intersections, inst->numIntersections, interSetStart, inst->interSetLen, numInterSet);
    
    // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fclose(f);
}



/**
 * Compute the frequencies of each variable inside the original constraints.
 * Saves array containing the variables indexes sorted by their frequencies in inst.
 *
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 *
 */
void computeVariableFrequency(int *izero, int *indexes, int nnz, instance *inst)
{
    int *frequencies = (int *)calloc(inst->num_cols, sizeof(int));
    int *variables = (int *)calloc(inst->num_cols, sizeof(int));

    // cycle the constraints, choosing the first one
    for (int i = 0; i < nnz; i++)
    {
        frequencies[indexes[i]]++;
    }

    // print_array_int(frequencies, inst->num_cols);
    int argMax;
    int Max;

    for (int i = 0; i < inst->num_cols; i++)
    {
        argMax = 0;
        Max = 0;

        for (int j = 0; j < inst->num_cols; j++)
        {
            if (frequencies[j] > Max)
            {
                Max = frequencies[j];
                argMax = j;
            }
        }
        variables[i] = argMax;
        frequencies[argMax] = 0;
    }
    inst->variableFreq = variables;
    free(frequencies);
}
