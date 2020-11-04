/**
 * This file contains functions needed for computing the solutions to the
 * scp using the CPLEX library.
 *
 * @author Petrucci Enrico
*/


#include "scp.h"


/** 
 * Callback function only useful to save the preprocessed redced problem as a new one.
 * 
 * @param type specifies the type of branch CPX_TYPE_VAR(0): variable branch. CPX_TYPE_SOS1(1): SOS1 branch. CPX_TYPE_SOS2(2): SOS2 branch. CPX_TYPE_ANY('A'): multiple bound changes or constraints will be used for branching
 * @param sos Specifies the special ordered set
 * @param nodecnt specifies the number og nodes cplex will create
 * @param bdcnt number of bounds changes defined in indices, lu, bd 
 * @param nodebeg discriminate in indices, lu, bd the changes in the different nodes
 * @param indices index of the variable
 * @param lu specifies if we add a lower or upper bound for the variable
 * @param bd specidies the new value of the bound
 * @param nodeest contains the estimations for the integer objective-function value the created nodes
 * @param useraction_p specifies if the branching was set by the user or if the one to use is the default one
 * 
 */
int legacyBranchingCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                         const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
{
   double start = second();
   instance *inst = (instance *)cbhandle; // casting of cbhandle
   int threadNum;
   CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &threadNum);
      
   //printf("callback\n");   
   if(inst->branching==0)
   {
      inst->defaultBranching[threadNum]++;
      double end = second();
      inst->timeInCallback[threadNum] += end - start;
      return 0;
   }
   else
   {
      if(nodecnt == 2)
      {
         int foundConstraint = 0;
         int i = -1;
         int j = -1;
         double *x = (double *) calloc(inst->num_cols, sizeof(double));
         double obj;

         // int node;
         // CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
            
         double* pseudocostDown = (double *) calloc(inst->num_cols, sizeof(double));
         double* pseudocostUp = (double *) calloc(inst->num_cols, sizeof(double));
         CPXgetcallbackpseudocosts(env, cbdata, wherefrom, pseudocostUp, pseudocostDown, 0, inst->num_cols-1);
            
         CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);
         //printf("LP relaxation solved optimally it has the objective %f\n", obj);
         CPXgetcallbacknodex(env, cbdata, wherefrom, x,  0, inst->num_cols-1);
         double Max = 0;
         double sec1 = second();
         int reducedConstraint = 1;

         if(reducedConstraint)
         {
            
            Max = findBestBranchingConstraintContainingVar(&i, indices[0], x, pseudocostDown, pseudocostUp, inst);
            inst->timeFindingConstraint[threadNum] += second()-sec1;
            // printf("comparing max constraint = %f and max variable = %f\n", Max, pseudocostDown[indices[0]]*x[indices[0]]*pseudocostUp[indices[0]]*(1-x[indices[0]]));
            // printf("First child:\n");
            // print_array_int(indices, nodebeg[1]);
            // print_array(bd, nodebeg[1]);
            // print_array_char(lu, nodebeg[1]);

            // printf("Second child:\n");
            // print_array_int(indices, bdcnt);
            // print_array(bd, bdcnt);
            // print_array_char(lu, bdcnt);
            double varCostDown = (0.001 > pseudocostDown[indices[0]]*x[indices[0]] ? 0.001 : pseudocostDown[indices[0]]*x[indices[0]]);
            double varCostUp = (0.001 > pseudocostUp[indices[0]]*(1-x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]]*(1-x[indices[0]]));
            
            if (i!= -1 && Max > 5*(varCostDown*varCostUp))
            {
                addBrachingChildsReduced(env, cbdata, wherefrom, obj, indices[0], i, inst);
                *useraction_p = CPX_CALLBACK_SET;
                inst->constraintBranching[threadNum]++;
                double end = second();
                inst->timeInCallback[threadNum] += end - start;
                return 0;
            }
            else
            {
                inst->defaultBranching[threadNum]++;
                double end = second();
                inst->timeInCallback[threadNum] += end - start;
                return 0;
            }
         }
         else
         {
            Max = findBestBranchingConstraint(&i, &j, x, pseudocostDown, pseudocostUp, inst);
            inst->timeFindingConstraint[threadNum] += second()-sec1;
            // printf("comparing max constraint = %f and max variable = %f\n", Max, pseudocostDown[indices[0]]*x[indices[0]]*pseudocostUp[indices[0]]*(1-x[indices[0]]));
            // printf("First child:\n");
            // print_array_int(indices, nodebeg[1]);
            // print_array(bd, nodebeg[1]);
            // print_array_char(lu, nodebeg[1]);

            // printf("Second child:\n");
            // print_array_int(indices, bdcnt);
            // print_array(bd, bdcnt);
            // print_array_char(lu, bdcnt);
            double varCostDown = (0.001 > pseudocostDown[indices[0]]*x[indices[0]] ? 0.001 : pseudocostDown[indices[0]]*x[indices[0]]);
            double varCostUp = (0.001 > pseudocostUp[indices[0]]*(1-x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]]*(1-x[indices[0]]));
            
            if (i!= -1 && Max > 2*(varCostDown*varCostUp))
            {
                addBrachingChilds(env, cbdata, wherefrom, obj, i, j, inst);
                *useraction_p = CPX_CALLBACK_SET;
                inst->constraintBranching[threadNum]++;
                double end = second();
                inst->timeInCallback[threadNum] += end - start;
                return 0;
            }
            else
            {
                inst->defaultBranching[threadNum]++;
                double end = second();
                inst->timeInCallback[threadNum] += end - start;
                return 0;
            }
         }
      }
      else
      {
         inst->defaultBranching[threadNum]++;
         double end = second();
         inst->timeInCallback[threadNum] += end - start;
         return 0;
      }
   }
}





/**
 * Solves an instance of the scp using the CPLEX library.
 *
 *
 * It creates the CPLEX environment and instantiate the CPLEX problem.
 * It then builds the model, solves the problem, gets the cost
 * of the optimal solution found, retrieves and saves this
 * information in the problem instance and finally frees and closes
 * the CPLEX model.
 *
 * @param inst instance of the scp
 * @returns 0 if executed successfully
*/
int scpopt(instance *inst)
{

   int error; // error = 0 -> no error

   /*
   * call to the CPLEX Callable Library, returns a pointer to the CPLEX environment
   * contains execution parameters (eg. timelimit, verbosity, etc.)
   */
   CPXENVptr env = CPXopenCPLEX(&error);

   /*
   * instantiates the problem objects, returns a pointer to the problem object
   * contains problem data (included the solution) there can be many problem objects
   */
   CPXLPptr lp = CPXcreateprob(env, &error, "scp");

   //CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL); // disable dynamic search

   // set random seed
   CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);

   // set time limit based on value passed as command line argument
   CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

   // import problem from file
   CPXreadcopyprob(env, lp, inst->input_file, "LP");

   inst->num_cols = CPXgetnumcols(env, lp);
   inst->num_rows = CPXgetnumrows(env, lp);

   char str[80];

   if(inst->callback==0)
   {
      solveUsingGenericCallback(env, lp, inst);
   }
   else
   {
      solveUsingLegacyCallback(env, lp, inst);
   }
   // free and close CPLEX model
   CPXfreeprob(env, &lp);
   CPXcloseCPLEX(&env);
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
void populateIntesectionsOf2(int* izero, int* indexes, int nnz, instance* inst)
{
   inst->intersections = (int **)calloc(inst->num_rows*(inst->num_rows-1)/2, sizeof(int*));
   inst->interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
   inst->interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
   inst->numInterSet = 0;
   inst->numIntersections = 0;

   // cycle the constraints, choosing the first one
   for(int i = 0; i < inst->num_rows-1; i++)
   {
      //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
      //Select the second constraint
      for(int n = i+1; n < inst->num_rows; n++)
      {
         int j = 0;
         int m = 0;
         int intersectionLen = 0;

         int length = (n<inst->num_rows-1 ? izero[n+1]-izero[n] : nnz-izero[n]);
         //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

         // alloc space for current intersection
         inst->intersections[inst->numIntersections] = (int *) calloc(inst->num_cols, sizeof(int));

         // cycle on all the variables inside both constraints
         while(j < izero[i+1]-izero[i] && m < length)// problem for the last constraint's length
         {
            // case 1: same variable. The variable is added to the computed information
            if(indexes[izero[i]+j]==indexes[izero[n]+m])
            {
               // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
               inst->intersections[inst->numIntersections][intersectionLen]=indexes[izero[i]+j];
               intersectionLen++;
               m++;
               j++;
            }
            // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
            if(indexes[izero[i]+j]<indexes[izero[n]+m])
            {
               j++;
            }

            if(indexes[izero[n]+m]<indexes[izero[i]+j])
            {
               m++;
            }
         }
         // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
         if(intersectionLen>1)
         {
            // realloc should free the unused spaced of the array
            inst->intersections[inst->numIntersections] = (int *) realloc(inst->intersections[inst->numIntersections], intersectionLen*sizeof(int));
            // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
            // case 1 validIntersection is empty
            if(inst->numIntersections==0)
            {
               // printf("Case 1 empty list\n");
               inst->interSetLen[0]=intersectionLen;
               inst->interSetStart[0]=0;
               inst->numInterSet++;
            }
            else
            {
               // cycle on all sets
               for(int k=0; k < inst->numInterSet; k++)
               {
                  // Case 2 Add new set at the beginning or add new set in between two sets
                  if((k==0 && intersectionLen > inst->interSetLen[k]) || (k>0 && intersectionLen<inst->interSetLen[k-1] && intersectionLen>inst->interSetLen[k]))
                  {
                     // printf("Case 2 new set in between\n");

                     // Update information on inst->interSetLen and inst->interSetStart
                     for(int a=inst->numInterSet; a>=k; a--)
                     {
                        inst->interSetLen[a+1]=inst->interSetLen[a];
                        inst->interSetStart[a+1]=inst->interSetStart[a]+1;
                     }

                     inst->interSetLen[k]=intersectionLen;
                     inst->numInterSet++;

                     // Updae infromation in validIntersection on the right with respect to k
                     for(int a=k; a<inst->numInterSet-1; a++)
                     {
                        // printf("from last position to %d\n", inst->interSetStart[a+1]-1);
                        int *temp=inst->intersections[inst->interSetStart[a+1]-1];
                        inst->intersections[inst->interSetStart[a+1]-1] = inst->intersections[inst->numIntersections];
                        inst->intersections[inst->numIntersections] = temp;
                     }
                     break;
                  }

                  // Case 3 no new set needed, only add an intersection and update everithing accordingly
                  if(intersectionLen == inst->interSetLen[k])
                  {

                     // printf("Case 3 same set\n");
                     for(int a=inst->numInterSet; a>k; a--)
                     {
                        inst->interSetStart[a]=inst->interSetStart[a]+1;
                     }

                     for(int a=k+1; a<inst->numInterSet; a++)
                     {
                        int *temp=inst->intersections[inst->interSetStart[a]-1];
                        // printf("inst->interSetStart[a]-1 = %d",inst->interSetStart[a]-1);
                        inst->intersections[inst->interSetStart[a]-1] = inst->intersections[inst->numIntersections];
                        inst->intersections[inst->numIntersections] = temp;
                     }
                     break;
                  }

                  // Case 4 Add new set at the end
                  if(k==inst->numInterSet-1 && intersectionLen < inst->interSetLen[k])
                  {

                     //printf("Case 4 new set at the end\n");
                     inst->interSetLen[inst->numInterSet] = intersectionLen;
                     inst->interSetStart[inst->numInterSet] = inst->numIntersections;
                     inst->numInterSet++;
                     break;
                  }
               }
            }
            inst->numIntersections++;
            // printf("Updated len and start:\n");
            // print_array_int(inst->interSetLen, inst->numInterSet);
            // print_array_int(inst->interSetStart, inst->numInterSet);
            // print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
            // printf("Total length %d\n", inst->numIntersections);
            // printf("breakpoint\n");
         }
         else
         {  // if the intersection is empty: free allocated space
            //printf("free\n");
            free(inst->intersections[inst->numIntersections]);
         }
      }
   }
   printf("Length for each set");
   print_array_int(inst->interSetLen, inst->numInterSet);
   printf("Starting index for each set for each set");
   print_array_int(inst->interSetStart, inst->numInterSet);
   //print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);

   printf("Total length %d max length %d\n", inst->numIntersections, inst->num_rows*(inst->num_rows-1)/2);

   FILE *f;
   f = fopen("Intersections.txt", "w");

   fprint_array_int_int(f, inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
   fclose(f);
   inst->shortestConstraint = inst->interSetLen[0];
   inst->lowestNumVariables = inst->interSetLen[0];
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
void populateIntesectionsOf2NoDup(int* izero, int* indexes, int nnz, instance* inst)
{
   inst->intersections = (int **)calloc(inst->num_rows*(inst->num_rows-1)/2, sizeof(int*));
   inst->interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
   inst->interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
   inst->numInterSet = 0;
   inst->numIntersections = 0;

   // cycle the constraints, choosing the first one
   for(int i = 0; i < inst->num_rows-1; i++)
   {
      //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
      //Select the second constraint
      for(int n = i+1; n < inst->num_rows; n++)
      {
         int j = 0;
         int m = 0;
         int intersectionLen = 0;

         int length = (n<inst->num_rows-1 ? izero[n+1]-izero[n] : nnz-izero[n]);
         //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

         // alloc space for current intersection
         inst->intersections[inst->numIntersections] = (int *) calloc(inst->num_cols, sizeof(int));

         // cycle on all the variables inside both constraints
         while(j < izero[i+1]-izero[i] && m < length)// problem for the last constraint's length
         {
            // case 1: same variable. The variable is added to the computed information
            if(indexes[izero[i]+j]==indexes[izero[n]+m])
            {
               // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
               inst->intersections[inst->numIntersections][intersectionLen]=indexes[izero[i]+j];
               intersectionLen++;
               m++;
               j++;
            }
            // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
            if(indexes[izero[i]+j]<indexes[izero[n]+m])
            {
               j++;
            }

            if(indexes[izero[n]+m]<indexes[izero[i]+j])
            {
               m++;
            }
         }
         // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
         if(intersectionLen>1)
         {
            // realloc should free the unused spaced of the array
            inst->intersections[inst->numIntersections] = (int *) realloc(inst->intersections[inst->numIntersections], intersectionLen*sizeof(int));
            // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
            // case 1 validIntersection is empty
            if(inst->numIntersections==0)
            {
               // printf("Case 1 empty list\n");
               inst->interSetLen[0]=intersectionLen;
               inst->interSetStart[0]=0;
               inst->numInterSet++;
               inst->numIntersections++;
            }
            else
            {
               //check if it is already present
               int alreadyIn = 0;
               for(int k=0; k < inst->numInterSet; k++)
               {
                  // if a subset has the same length
                  if(inst->interSetLen[k]==intersectionLen)
                  {
                     // for each intersection in the subset
                     for(int a = inst->interSetStart[k]; a < (k == inst->numInterSet-1 ? inst->numIntersections : inst->interSetStart[k+1]); a++ )
                     {
                        int counter = 0;
                        //for each variable in a single subset
                        for(int b = 0; b < inst->interSetLen[k]; b++)
                        {
                           if(inst->intersections[a][b] != inst->intersections[inst->numIntersections][b])
                           {
                              break;
                           }
                           else
                           {
                              counter++;
                           }
                        }
                        if(counter==inst->interSetLen[k])
                        {
                           alreadyIn=1;
                           break;
                        }
                     }
                     break;
                  }
               }
               if(alreadyIn==0)
               {
                  // cycle on all sets
                  for(int k=0; k < inst->numInterSet; k++)
                  {
                     // Case 2 Add new set at the beginning or add new set in between two sets
                     if((k==0 && intersectionLen > inst->interSetLen[k]) || (k>0 && intersectionLen<inst->interSetLen[k-1] && intersectionLen>inst->interSetLen[k]))
                     {
                        // printf("Case 2 new set in between\n");

                        // Update information on inst->interSetLen and inst->interSetStart
                        for(int a=inst->numInterSet; a>=k; a--)
                        {
                           inst->interSetLen[a+1]=inst->interSetLen[a];
                           inst->interSetStart[a+1]=inst->interSetStart[a]+1;
                        }

                        inst->interSetLen[k]=intersectionLen;
                        inst->numInterSet++;

                        // Updae infromation in validIntersection on the right with respect to k
                        for(int a=k; a<inst->numInterSet-1; a++)
                        {
                           // printf("from last position to %d\n", inst->interSetStart[a+1]-1);
                           int *temp=inst->intersections[inst->interSetStart[a+1]-1];
                           inst->intersections[inst->interSetStart[a+1]-1] = inst->intersections[inst->numIntersections];
                           inst->intersections[inst->numIntersections] = temp;
                        }
                        break;
                     }

                     // Case 3 no new set needed, only add an intersection and update everithing accordingly
                     if(intersectionLen == inst->interSetLen[k])
                     {

                        // printf("Case 3 same set\n");
                        for(int a=inst->numInterSet; a>k; a--)
                        {
                           inst->interSetStart[a]=inst->interSetStart[a]+1;
                        }

                        for(int a=k+1; a<inst->numInterSet; a++)
                        {
                           int *temp=inst->intersections[inst->interSetStart[a]-1];
                           // printf("inst->interSetStart[a]-1 = %d",inst->interSetStart[a]-1);
                           inst->intersections[inst->interSetStart[a]-1] = inst->intersections[inst->numIntersections];
                           inst->intersections[inst->numIntersections] = temp;
                        }
                        break;
                     }

                     // Case 4 Add new set at the end
                     if(k==inst->numInterSet-1 && intersectionLen < inst->interSetLen[k])
                     {

                        //printf("Case 4 new set at the end\n");
                        inst->interSetLen[inst->numInterSet] = intersectionLen;
                        inst->interSetStart[inst->numInterSet] = inst->numIntersections;
                        inst->numInterSet++;
                        break;
                     }
                  }
                  inst->numIntersections++;
               }
            }
            
            // printf("Updated len and start:\n");
            // print_array_int(inst->interSetLen, inst->numInterSet);
            // print_array_int(inst->interSetStart, inst->numInterSet);
            // print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
            // printf("Total length %d\n", inst->numIntersections);
            // printf("breakpoint\n");
         }
         else
         {  // if the intersection is empty: free allocated space
            //printf("free\n");
            free(inst->intersections[inst->numIntersections]);
         }
      }
   }
   printf("Length for each set");
   print_array_int(inst->interSetLen, inst->numInterSet);
   printf("Starting index for each set for each set");
   print_array_int(inst->interSetStart, inst->numInterSet);
   //print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);

   printf("Total length %d max length %d\n", inst->numIntersections, inst->num_rows*(inst->num_rows-1)/2);

   FILE *f;
   f = fopen("IntersectionsNoDup.txt", "w");

   fprint_array_int_int(f, inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
   fclose(f);
   inst->shortestConstraint = inst->interSetLen[0];
   inst->lowestNumVariables = inst->interSetLen[0];
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
void populateIntesectionsOf3(int* izero, int* indexes, int nnz, instance* inst)
{
   int possibleIntersections = (inst->num_rows)*(inst->num_rows-1)*(inst->num_rows-2)/6;
   inst->intersections = (int **)calloc(possibleIntersections, sizeof(int*));
   inst->interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
   inst->interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
   inst->numInterSet = 0;
   inst->numIntersections = 0;

   // cycle the constraints, choosing the first one
   for(int x = 0; x < inst->num_rows-2; x++)
   {
      //printf("x = %d\n", x);
      for(int i = x+1; i < inst->num_rows-1; i++)
      {
         //printf("i = %d\n",i);
         //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
         //Select the second constraint
         for(int n = i+1; n < inst->num_rows; n++)
         {
            //printf("n = %d\n",n);
            int y = 0;
            int j = 0;
            int m = 0;
            int intersectionLen = 0;

            int length = (n<inst->num_rows-1 ? izero[n+1]-izero[n] : nnz-izero[n]);
            //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);

            // alloc space for current intersection
            inst->intersections[inst->numIntersections] = (int *) calloc(inst->num_cols, sizeof(int));

            // cycle on all the variables inside both constraints
            while(y < izero[x+1]-izero[x] && j < izero[i+1]-izero[i] && m < length)// problem for the last constraint's length
            {
               //printf("y = %d, j = %d, m = %d\n", y, j, m);
               //printf("first variable = %d, second variable = %d, third variable = %d\n", indexes[izero[x]+y], indexes[izero[i]+j], indexes[izero[n]+m]);
               // case 1: same variable. The variable is added to the computed information
               if(indexes[izero[x]+y] == indexes[izero[i]+j] && indexes[izero[i]+j]==indexes[izero[n]+m])
               {
                  //printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                  inst->intersections[inst->numIntersections][intersectionLen]=indexes[izero[i]+j];
                  intersectionLen++;
                  y++;
                  m++;
                  j++;

               }
               if(indexes[izero[x]+y]<indexes[izero[i]+j] && indexes[izero[x]+y]<indexes[izero[n]+m])
               {
                  y++;
               }
               if(indexes[izero[i]+j]<indexes[izero[x]+y] && indexes[izero[i]+j]<indexes[izero[n]+m])
               {
                  j++;
               }
               if(indexes[izero[n]+m]<indexes[izero[x]+y] && indexes[izero[n]+m]<indexes[izero[i]+j])
               {
                  m++;
               }
               if(indexes[izero[x]+y]>indexes[izero[i]+j] && indexes[izero[x]+y]>indexes[izero[n]+m])
               {
                  j++;
                  m++;
               }
               if(indexes[izero[i]+j]>indexes[izero[x]+y] && indexes[izero[i]+j]>indexes[izero[n]+m])
               {
                  y++;
                  m++;
               }
               if(indexes[izero[n]+m]>indexes[izero[x]+y] && indexes[izero[n]+m]>indexes[izero[i]+j])
               {
                  y++;
                  j++;
               }
            }
            // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
            if(intersectionLen>1)
            {
               // realloc should free the unused spaced of the array
               inst->intersections[inst->numIntersections] = (int *) realloc(inst->intersections[inst->numIntersections], intersectionLen*sizeof(int));
               //printf("Intersection between constraints %d, %d and %d has %d variables\n", x, i, n, intersectionLen);
               // case 1 validIntersection is empty
               if(inst->numIntersections==0)
               {
                  // printf("Case 1 empty list\n");
                  inst->interSetLen[0]=intersectionLen;
                  inst->interSetStart[0]=0;
                  inst->numInterSet++;
               }
               else
               {
                  // cycle on all sets
                  for(int k=0; k < inst->numInterSet; k++)
                  {
                     // Case 2 Add new set at the beginning or add new set in between two sets
                     if((k==0 && intersectionLen > inst->interSetLen[k]) || (k>0 && intersectionLen<inst->interSetLen[k-1] && intersectionLen>inst->interSetLen[k]))
                     {
                        // printf("Case 2 new set in between\n");

                        // Update information on inst->interSetLen and inst->interSetStart
                        for(int a=inst->numInterSet; a>=k; a--)
                        {
                           inst->interSetLen[a+1]=inst->interSetLen[a];
                           inst->interSetStart[a+1]=inst->interSetStart[a]+1;
                        }

                        inst->interSetLen[k]=intersectionLen;
                        inst->numInterSet++;

                        // Update infromation in validIntersection on the right with respect to k
                        for(int a=k; a<inst->numInterSet-1; a++)
                        {
                           // printf("from last position to %d\n", inst->interSetStart[a+1]-1);
                           int *temp=inst->intersections[inst->interSetStart[a+1]-1];
                           inst->intersections[inst->interSetStart[a+1]-1] = inst->intersections[inst->numIntersections];
                           inst->intersections[inst->numIntersections] = temp;
                        }
                        break;
                     }

                     // Case 3 no new set needed, only add an intersection and update everithing accordingly
                     if(intersectionLen == inst->interSetLen[k])
                     {

                        // printf("Case 3 same set\n");
                        for(int a=inst->numInterSet; a>k; a--)
                        {
                           inst->interSetStart[a]=inst->interSetStart[a]+1;
                        }

                        for(int a=k+1; a<inst->numInterSet; a++)
                        {
                           int *temp=inst->intersections[inst->interSetStart[a]-1];
                           // printf("inst->interSetStart[a]-1 = %d",inst->interSetStart[a]-1);
                           inst->intersections[inst->interSetStart[a]-1] = inst->intersections[inst->numIntersections];
                           inst->intersections[inst->numIntersections] = temp;
                        }
                        break;
                     }

                     // Case 4 Add new set at the end
                     if(k==inst->numInterSet-1 && intersectionLen < inst->interSetLen[k])
                     {

                        //printf("Case 4 new set at the end\n");
                        inst->interSetLen[inst->numInterSet] = intersectionLen;
                        inst->interSetStart[inst->numInterSet] = inst->numIntersections;
                        inst->numInterSet++;
                        break;
                     }
                  }
               }
               inst->numIntersections++;
               // printf("Updated len and start:\n");
               // print_array_int(inst->interSetLen, inst->numInterSet);
               // print_array_int(inst->interSetStart, inst->numInterSet);
               // print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
               // printf("Total length %d\n", inst->numIntersections);
               // printf("breakpoint\n");
            }
            else
            {  // if the intersection is empty: free allocated space
               //printf("free\n");
               free(inst->intersections[inst->numIntersections]);
            }
         }
      }
   }
   printf("Length for each set");
   print_array_int(inst->interSetLen, inst->numInterSet);
   printf("Starting index for each set for each set");
   print_array_int(inst->interSetStart, inst->numInterSet);
   //print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);

   printf("Total length %d max length %d\n", inst->numIntersections, possibleIntersections);

   //FILE *f;
   //f = fopen("Intersections.txt", "w");

   //fprint_array_int_int(f, inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
   //fclose(f);
   inst->shortestConstraint = inst->interSetLen[0];
   inst->lowestNumVariables = inst->interSetLen[0];
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
void computeVariableFrequency(int* izero, int* indexes, int nnz, instance* inst)
{
   int* frequencies = (int *)calloc(inst->num_cols, sizeof(int));
   int* variables = (int *)calloc(inst->num_cols, sizeof(int));

   // cycle the constraints, choosing the first one
   for(int i = 0; i < nnz; i++)
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

      for(int j = 0; j < inst->num_cols; j++)
      {
         if(frequencies[j]>Max)
         {
            Max=frequencies[j];
            argMax=j;
         }
      }
      variables[i]=argMax;
      frequencies[argMax]=0;
   }
   inst->variableFreq = variables;
   free(frequencies);
}


double findBestBranchingConstraint(int *n, int *m, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst)
{
   double* estimateScoreDown = (double *) calloc(inst->numIntersections, sizeof(double));
   double* estimateScoreUp = (double *) calloc(inst->numIntersections, sizeof(double));
   double* constraintScorePrevision = (double *) calloc(inst->numIntersections, sizeof(double));
   double epsilon = 0.0001;

   double Max = 0;
   int best_i=-1;
   int best_j=-1;

    
    for(int i=0; i < inst->numInterSet; i++)
    {
        double pseudoDown;
        double pseudoUp;
        // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
        // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
        int j=0;
        // cycle on all the intersection having the same lenght
        for (; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
        {
            // printf("j=%d\n", j);
            double sum = 0;
            // cycle all the variables in the intersection and get the sum
            pseudoDown=0;
            pseudoUp=INT_MAX; 
            
            // printf("***Managing new constraint %d ***\n", inst->interSetStart[i]+j);
            for(int k = 0; k < inst->interSetLen[i]; k++)
            {
                // printf("variable %d has value %f and down pseudocost = %f\n", inst->intersections[inst->interSetStart[i]+j][k], rootSolution[inst->intersections[inst->interSetStart[i]+j][k]], pseudocostDown[inst->intersections[inst->interSetStart[i]+j][k]]);
                //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
                // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
                sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
                pseudoDown += pseudocostDown[inst->intersections[inst->interSetStart[i]+j][k]]*x[inst->intersections[inst->interSetStart[i]+j][k]];
                // printf("current downpseudocost = %f\n", pseudoDown);
                double candidate = pseudocostUp[inst->intersections[inst->interSetStart[i]+j][k]];
                
                pseudoUp = (pseudoUp > candidate ? candidate : pseudoUp);
            }
            
            pseudoUp *= (1-sum);

            if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i]) // set covering
            {
                // printf("breakpoint");
                estimateScoreDown[inst->interSetStart[i]+j] = pseudoDown;
                estimateScoreUp[inst->interSetStart[i]+j] = pseudoUp;
                constraintScorePrevision[inst->interSetStart[i]+j] = pseudoDown*pseudoUp;
            }
            else
            {
                estimateScoreDown[inst->interSetStart[i]+j] = epsilon;
                estimateScoreUp[inst->interSetStart[i]+j] = epsilon;
                
                //printf("Constraint %d  not usable\n",inst->interSetStart[i]+j);
                constraintScorePrevision[inst->interSetStart[i]+j] = epsilon*epsilon;
            }
            //print_array(productScoreConstraint, inst->interSetStart[i]+j+1);
            //printf("breakpoint\n");
            //printf("score: %f\n", productScoreConstraint[inst->interSetStart[i]+j]);
            if(Max<constraintScorePrevision[inst->interSetStart[i]+j])
            {
                Max=constraintScorePrevision[inst->interSetStart[i]+j];
                best_i=i;
                best_j=j;
                //printf("new max for Constraints prevision: %f\n", Max);
            }
        }
    }
    *n = best_i;
    *m = best_j;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return(Max);
}
    

double findBestBranchingConstraintContainingVar(int *n, int bestVar, double* x, double* pseudocostDown, double* pseudocostUp, instance* inst)
{

   double* estimateScoreDown = (double *) calloc(inst->constraintCounter[bestVar], sizeof(double));
   double* estimateScoreUp = (double *) calloc(inst->constraintCounter[bestVar], sizeof(double));
   double* constraintScorePrevision = (double *) calloc(inst->constraintCounter[bestVar], sizeof(double));
   double epsilon = 0.0001;

   double Max = 0;
   int best_i = -1;
    
    for(int i=0; i < inst->constraintCounter[bestVar]; i++)
    {
        double pseudoDown;
        double pseudoUp;
        
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        pseudoDown=0;
        pseudoUp=INT_MAX; 
        
        // printf("***Managing new constraint %d ***\n", inst->interSetStart[i]+j);
        for(int k = 0; k < inst->varConstrDim[bestVar][i]; k++)
        {
            // printf("variable %d has value %f and down pseudocost = %f\n", inst->intersections[inst->interSetStart[i]+j][k], rootSolution[inst->intersections[inst->interSetStart[i]+j][k]], pseudocostDown[inst->intersections[inst->interSetStart[i]+j][k]]);
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            int* constraint = inst->intersections[inst->varConstrTable[bestVar][i]];
            sum+=x[constraint[k]];
            pseudoDown += pseudocostDown[constraint[k]]*x[constraint[k]];
            // printf("current downpseudocost = %f\n", pseudoDown);
            double candidate = pseudocostUp[constraint[k]];
            
            pseudoUp = (pseudoUp > candidate ? candidate : pseudoUp);
        }
        
        pseudoUp *= (1-sum);

        if(sum>1e-05*inst->varConstrDim[bestVar][i] && sum<1-1e-05*inst->varConstrDim[bestVar][i]) // set covering
        {
            // printf("breakpoint");
            estimateScoreDown[i] = pseudoDown;
            estimateScoreUp[i] = pseudoUp;
            constraintScorePrevision[i] = pseudoDown*pseudoUp;
        }
        else
        {
            estimateScoreDown[i] = epsilon;
            estimateScoreUp[i] = epsilon;
            
            //printf("Constraint %d  not usable\n",inst->interSetStart[i]+j);
            constraintScorePrevision[i] = epsilon*epsilon;
        }
        //print_array(productScoreConstraint, inst->interSetStart[i]+j+1);
        //printf("breakpoint\n");
        //printf("score: %f\n", productScoreConstraint[inst->interSetStart[i]+j]);
        if(Max<constraintScorePrevision[i])
        {
            Max=constraintScorePrevision[i];
            best_i=i;
            //printf("new max for Constraints prevision: %f\n", Max);
        }
    }
    *n = best_i;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return(Max);
}



void addBrachingChilds(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, int j, instance* inst)
{
         
   // generating information for adding the constraint to CPLEX
   // in the first case the sum of the variables is set equal to zero,
   // in the second case the sum of the variables is set to be equal or greater than 1
   double rhs1 = 0;
   double rhs2 = 1;
   char sense1 = 'E';
   char sense2 = 'G';
   int izero = 0;
   // generating the array of values, all ones for the variables selected
   double* values = (double *) calloc(inst->interSetLen[i], sizeof(double));
   for(int n=0; n<inst->interSetLen[i]; n++)
   {
      values[n]=1;
   }

   /* We want to branch. */
   int child1, child2;
   // print_array_int(inst->intersections[inst->interSetStart[i]+j], inst->interSetLen[i]);
   // print_array(values, inst->interSetLen[i]);
   // printf("adding constraints with %d variables, rhs1=%f, rhs2=%f, sense1=%c, sense2=%c, izero=%d\n", inst->interSetLen[i], rhs1, rhs2, sense1, sense2, izero);
   if(inst->interSetLen[i]<inst->shortestConstraint)
   {
      inst->shortestConstraint = inst->interSetLen[i];
      printf("branching on a constraint with %d variables\n", inst->interSetLen[i]);

   }
   if(CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->interSetLen[i], &rhs1, "E", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, NULL, &child1)!=0)
   {
      print_error("First branch was not succesful\n");
   }
   //if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs2, "E", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child2)==0)
   if(CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->interSetLen[i], &rhs2, "G", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, NULL, &child2)!=0)
   {
      print_error("Second branch was not succesful\n");
   }
   free(values);
}


void addBrachingChildsReduced(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int bestVar, int i, instance* inst)
{
         
   // generating information for adding the constraint to CPLEX
   // in the first case the sum of the variables is set equal to zero,
   // in the second case the sum of the variables is set to be equal or greater than 1
   double rhs1 = 0;
   double rhs2 = 1;
   char sense1 = 'E';
   char sense2 = 'G';
   int izero = 0;
   // generating the array of values, all ones for the variables selected
   double* values = (double *) calloc(inst->varConstrDim[bestVar][i], sizeof(double));
   for(int n=0; n<inst->varConstrDim[bestVar][i]; n++)
   {
      values[n]=1;
   }

   /* We want to branch. */
   int child1, child2;
   // print_array_int(inst->intersections[inst->interSetStart[i]+j], inst->interSetLen[i]);
   // print_array(values, inst->interSetLen[i]);
   // printf("adding constraints with %d variables, rhs1=%f, rhs2=%f, sense1=%c, sense2=%c, izero=%d\n", inst->interSetLen[i], rhs1, rhs2, sense1, sense2, izero);
   
   if(CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->varConstrDim[bestVar][i], &rhs1, "E", &izero, inst->intersections[inst->varConstrTable[bestVar][i]], values, obj, NULL, &child1)!=0)
   {
      print_error("First branch was not succesful\n");
   }
   //if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs2, "E", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child2)==0)
   if(CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->varConstrDim[bestVar][i], &rhs2, "G", &izero, inst->intersections[inst->varConstrTable[bestVar][i]], values, obj, NULL, &child2)!=0)
   {
      print_error("Second branch was not succesful\n");
   }
   free(values);
}


void solveUsingLegacyCallback(CPXENVptr env, CPXLPptr lp, instance* inst)
{
   char str[80];

   if(inst->branching==0)
      sprintf(str, "logfile_defaultBranchinglegacy.txt");
   else
      sprintf(str, "logfile_constraintBranchinglegacy.txt");
   CPXsetlogfilename(env, str, "w");

   CPXsetbranchcallbackfunc(env, legacyBranchingCallback, inst);

   if (inst->threads==0)
      CPXgetnumcores(env, &inst->threads);
   CPXsetintparam(env, CPX_PARAM_THREADS, inst->threads);

   CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
   // CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);//CPX_VARSEL_PSEUDO);//
   CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_PSEUDO);

   inst->constraintBranching = (int *) calloc(inst->threads, sizeof(int));
   inst->defaultBranching = (int *) calloc(inst->threads, sizeof(int));

   inst->timeInCallback = (double *) calloc(inst->threads, sizeof(double));

   inst->timeFindingConstraint = (double *) calloc(inst->threads, sizeof(double));

   if(inst->branching == 1)
   {
      //CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);//CPX_VARSEL_PSEUDO);//
   
      // Get the rows of the problem
      int nnz;
      int *izero = (int *) calloc(inst->num_rows, sizeof(int));
      int *indexes = (int *) calloc(inst->num_rows*inst->num_cols, sizeof(int));// this one is big!
      double *values = (double *) calloc(inst->num_rows*inst->num_cols, sizeof(double));// this one is big too!
      int surplus_p;

      CPXgetrows(env, lp, &nnz, izero, indexes, values, inst->num_rows*inst->num_cols, &surplus_p, 0, inst->num_rows-1);

      indexes = (int *) realloc(indexes, nnz*sizeof(int));// realloc should free the unused spaced of the array
      values = (double *) realloc(values, nnz*sizeof(double));// realloc should free the unused spaced of the array
      double start = second();
      populateIntesectionsOf2(izero, indexes, nnz, inst);
      double end = second();

      printf("Found %d constraint intersections in %f\n", inst->numIntersections, end-start);
      
      start = second();
      populateVariableConstraintTable(inst);
      end = second();
      
      printf("Populated variable-constraint table in %f\n", end-start);
      
      

      // start = second();
      // populateIntesectionsOf3(izero, indexes, nnz, inst);
      // //populateIntesectionsOf4(inst);// not working yet
      // end = second();
      // printf("Found %d constraint intersections in %f\n", inst->numIntersections, end-start);


      start = second();
      computeVariableFrequency(izero, indexes, nnz, inst);
      //populateIntesectionsOf4(inst);// not working yet
      end = second();
      printf("Computed variable frequencies in %f\n", end-start);

      // print_array_int(inst->variableFreq, inst->num_cols);
      // int ncores = 1;
      // //CPXgetnumcores(env, &ncores);
      // CPXsetintparam(env, CPX_PARAM_THREADS, ncores); // it was reset after callback

      free(izero);
      free(indexes);
      free(values);
   }
   if (CPXmipopt(env, lp)) // solve the problem
   print_error(" Problems on CPXmipopt");

   inst->solution = (double *) calloc(inst->num_cols, sizeof(double));
   double obj;

   if(CPXgetx(env, lp, inst->solution, 0, inst->num_cols-1)==0 && CPXgetbestobjval(env, lp, &obj)==0)
      printf("Solution found with objval = %f\n", obj);
   else
      print_error("Problems in getting the solution");
   
   int constraintCount=0;
   int defaultCount=0;
   double callbackTime = 0;
   double findingConstraintTime = 0;

   if(inst->branching == 1)
   {
         for(int i=0; i<inst->numIntersections; i++)
         {
            free(inst->intersections[i]);
         }
         free(inst->intersections);
         free(inst->interSetLen);
         free(inst->interSetStart);
         free(inst->variableFreq);
   }
   for(int i=0; i<inst->threads; i++)
   {
      constraintCount += (inst->constraintBranching[i]);
      defaultCount += (inst->defaultBranching[i]);
      callbackTime += (inst->timeInCallback[i]);
      findingConstraintTime += (inst->timeFindingConstraint[i]);
   }
   printf("%d constraint branching\n%d default branching\n", constraintCount, defaultCount);
   
   
   printf("Total time in callback = %f which is %f%% of the total\n", callbackTime, 100*callbackTime/(second()-inst->startTime));
   printf("Time spent in callback choosing the best constraint: %f\n", findingConstraintTime);
   
   free(inst->solution);
   free(inst->constraintBranching);
   free(inst->defaultBranching);
}
        

void populateVariableConstraintTable(instance *inst)
{
    int ** varConstrTable = (int **) calloc(inst->num_cols, sizeof(int*));
    int* constraintCounter = (int *) calloc(inst->num_cols, sizeof(int)); 
    
    int ** varConstrDim = (int **) calloc(inst->num_cols, sizeof(int*));
 

    for(int i = 0; i<inst->num_cols; i++)
    {
        varConstrTable[i] = (int *) calloc(inst->numIntersections, sizeof(int));
        varConstrDim[i] = (int *) calloc(inst->numIntersections, sizeof(int));
    }
   
   // cycle on the set having the same length
    for(int i=0; i < inst->numInterSet; i++)
    {
        int j=0;
        // cycle on the constraint having the same length
        for (; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
        {
            // printf("analyzing constraint:");
            // print_array_int(inst->intersections[inst->interSetStart[i]+j], inst->interSetLen[i]);
            // cycle on the variables
            for(int k = 0; k < inst->interSetLen[i]; k++)
            {
                int currentVar = (inst->intersections[inst->interSetStart[i]+j][k]);
                varConstrTable[currentVar][constraintCounter[currentVar]] = inst->interSetStart[i]+j;
                varConstrDim[currentVar][constraintCounter[currentVar]] = inst->interSetLen[i];
                constraintCounter[currentVar]++;
            }
        }

    }
    FILE *f;
    f = fopen("varToConstraints.txt", "w");
    fprint_array_int_int1(f, varConstrTable, varConstrDim, constraintCounter, inst->num_cols);
    fclose(f);
    inst->varConstrTable = varConstrTable;
    inst->constraintCounter = constraintCounter; 
    inst->varConstrDim = varConstrDim;
}      