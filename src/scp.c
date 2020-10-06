/**
 * This file contains functions needed for computing the solutions to the
 * scp using the CPLEX library.
 * 
 * @author Petrucci Enrico
*/


#include "scp.h"

int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle) 
{
   if ( contextid == CPX_CALLBACKCONTEXT_BRANCHING)
   {
      instance *inst = (instance *)cbhandle; // casting of cbhandle
      printf("inside branching generic callback\n");
      /* Check the status of the continuous relaxation. */
      int statind;

      // with 0 it forces the optimality of the lp-relaxation
      CPXcallbackgetrelaxationstatus(context, &statind, 0);

      if (statind == CPX_STAT_OPTIMAL || statind == CPX_STAT_OPTIMAL_INFEAS)
      {
         double *x = (double *) calloc(inst->num_cols, sizeof(double));
         double obj;
         CPXcallbackgetrelaxationpoint(context, x, 0, inst->num_cols - 1, &obj);
         printf("LP relaxation solved optimally it has the objective %f\n", obj);
         int foundConstraint = 0;
         int i=0;
         int j=0;
         for(; i<inst->numInterSet; i++)
         {
            printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
            printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
            j=0;
            for (; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
            {
               printf("j=%d\n", j);
               double sum = 0;
               for(int k = 0; k < inst->interSetLen[i]; k++)
               {
                  //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
                  printf("number %d, variable %d has coefficient %f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
                  sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
               }
               printf("Constraint: i=%d, j=%d has sum in the solution = %f \n", i, j, sum);   
               if(sum!=0 && sum!=1)
               {
                  printf("Found constraint not yet used: i=%d, j=%d\n", i, j);
                  foundConstraint=1;
                  break;
               }
            }
            if(foundConstraint==1)
            {
               break;
            }
         }
         double rhs1 = 0;
         double rhs2 = 1;
         char sense = 'E';
         int izero = 0;
         double* values = (double *) calloc(inst->interSetLen[i], sizeof(double));
         for(int n=0; n<inst->interSetLen[i]; n++)
         {
            values[n]=1;
         }

         //    /* We want to branch. */
         int child1, child2;
         printf("i=%d, j=%d",i,j);
         print_array_int(inst->intersections[inst->interSetStart[i]+j], inst->interSetLen[i]);
         // if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs1, &sense, &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child1)==0)
         // {
         //    printf("First branch was succesful\n");
         // }
         // else
         // {
         //    printf("First branch was not succesful\n");
         // }
         // if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs2, &sense, &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child2)==0)
         // {
         //    printf("Second branch was succesful\n");
         // }
         // else
         // {
         //    printf("Second branch was not succesful\n");
         // }

         //free(values);
         CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, obj, &child1);
         CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 0, 0, NULL, NULL, NULL, NULL, NULL, obj, &child2);
         printf("breakpoint\n");
      }
      else 
      {
         printf("LP relaxation not solved optimally\n");
      }
   }
}

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

/**
 * Solves an instance of the scp using the CPLEX library.
 * 
 * 
 * It opens the CPLEX model, by getting a pointer to the CPLEX 
 * environment and instantiating the problem object.
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
   * call to the CPLEX Callable Library
   * returns a pointer to the CPLEX environment
   * contains execution parameters (eg. timelimit, verbosity, etc.)
   */
   CPXENVptr env = CPXopenCPLEX(&error);


   /*
   * instantiates the problem objects
   * returns a pointer to the problem object
   * contains problem data (included the solution) 
   * there can be many problem objects
   */

   CPXLPptr lp = CPXcreateprob(env, &error, "scp");

   //CPXsetintparam(env, CPX_PARAM_MIPSEARCH, CPX_MIPSEARCH_TRADITIONAL);
   // set random seed
   CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);

   // set time limit based on value passed as command line argument
   CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);

   CPXsetlogfilename(env, "logfile.txt", "w");

   CPXreadcopyprob(env, lp, inst->input_file, "LP");


   if(inst->extractPreprocessing==1)
   {
      // install legacy callback for extracting the preprocessed problem
      printf("Installing legacy callback for preprocessing node\n");
      //CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // let MIP callbacks work on the original model
      CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 1); // Set represolve always at the same level
      // disable all the cuts
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_BQP, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_Cliques, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_Covers, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_Disjunctive, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_FlowCovers, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_PathCut, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_Gomory, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_GUBCovers, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_Implied, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_LocalImplied, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_LiftProj, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_MIRCut, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_MCFCut	, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_RLT, -1);
      CPXsetintparam(env, CPXPARAM_MIP_Cuts_ZeroHalfCut, -1);
      


      CPXsetbranchcallbackfunc(env, preprocessinglegacycallback, inst );

      inst->num_cols = CPXgetnumcols(env, lp);
      inst->num_rows = CPXgetnumrows(env, lp);
      printf("The original models has %d columns\n", inst->num_cols);
      printf("The original models has %d rows\n", inst->num_rows);

      if (CPXmipopt(env, lp)) // solve the problem
         printf(" Problems on CPXmipopt");

      
   }
   else
   {
      if(inst->callback==0)
      {
         printf("Solving without using callbacks\n");
      }
      if(inst->callback==1)
      {
         inst->num_cols = CPXgetnumcols(env, lp);
         inst->num_rows = CPXgetnumrows(env, lp);
    
         int nnz;
         int *izero = (int *) calloc(inst->num_rows, sizeof(int));
         int *indexes = (int *) calloc(inst->num_rows*inst->num_cols, sizeof(int));// this one is big!
         double *values = (double *) calloc(inst->num_rows*inst->num_cols, sizeof(double));// this one is big too!
         int space;
         int surplus_p;


         CPXgetrows(env, lp, &nnz, izero, indexes, values, inst->num_rows*inst->num_cols, &surplus_p, 0, inst->num_rows-1);

         indexes = (int *) realloc(indexes, nnz*sizeof(int));// realloc should free the unused spaced of the array
         values = (double *) realloc(values, nnz*sizeof(double));// realloc should free the unused spaced of the array
         

         int ** validIntersections = (int **)calloc(inst->num_rows*(inst->num_rows-1)/2, sizeof(int*)); // superbig
         int * interSetLen = (int *)calloc(inst->num_cols, sizeof(int)); // superbig
         int * interSetStart = (int *)calloc(inst->num_cols, sizeof(int)); // superbig
         
         int numInterSet = 0;
         int numIntersection = 0;
         
         // select the first constraint
         for(int i = 0; i < inst->num_rows-1; i++)
         {
            int length = (i<inst->num_rows-1 ? (izero[i+1] - izero[i]) : (nnz - izero[i]));
            //printf("Considering i = %d, starts in %d it has length %d\n", i, izero[i], length);
            //Select the second constraint
            for(int n = i+1; n < inst->num_rows; n++)
            {
               int j = 0;
               int m = 0;
               int intersectionLen = 0;

               length = (n<inst->num_rows-1 ? izero[n+1]-izero[n] : nnz-izero[n]);
               //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);
               validIntersections[numIntersection] = (int *) calloc(inst->num_cols, sizeof(int));
               
               while(j < izero[i+1]-izero[i] && m < length)// problem for the last constraint's length
               {
                  if(indexes[izero[i]+j]==indexes[izero[n]+m])
                  {
                     // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
                  
                     validIntersections[numIntersection][intersectionLen]=indexes[izero[i]+j];
                     intersectionLen++;
                     m++;
                     j++;
                  }

                  if(indexes[izero[i]+j]<indexes[izero[n]+m])
                  {
                     j++;
                  }
                  
                  if(indexes[izero[n]+m]<indexes[izero[i]+j])
                  {
                     m++;
                  }
               }
               if(intersectionLen>0)
               {
                  validIntersections[numIntersection] = (int *) realloc(validIntersections[numIntersection], intersectionLen*sizeof(int));// realloc should free the unused spaced of the array
                  // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
                  if(numIntersection==0)
                  {
                     // printf("Case 1 empty list\n");
                     interSetLen[0]=intersectionLen;
                     interSetStart[0]=0;
                     numInterSet++;
                  }
                  else
                  {
                     for(int k=0; k < numInterSet; k++)
                     {
                        //Add new set at the beginning or add new set in between two sets
                        if((k==0 && intersectionLen > interSetLen[k]) || (k>0 && intersectionLen<interSetLen[k-1] && intersectionLen>interSetLen[k]))
                        {
                           // printf("Case 2 new set in between\n");
                           for(int a=numInterSet; a>=k; a--)
                           {
                              interSetLen[a+1]=interSetLen[a];
                              interSetStart[a+1]=interSetStart[a]+1;
                           }

                           interSetLen[k]=intersectionLen;
                           numInterSet++;
                
                           for(int a=k; a<numInterSet-1; a++)
                           {
                              // printf("from last position to %d\n", interSetStart[a+1]-1);
                              int *temp=validIntersections[interSetStart[a+1]-1];
                              validIntersections[interSetStart[a+1]-1] = validIntersections[numIntersection];
                              validIntersections[numIntersection] = temp;
                           }
                           //numInterSet++;
                           break;
                        }
                        
                        //No new set needed, only update
                        if(intersectionLen == interSetLen[k])
                        {

                           // printf("Case 3 same set\n");
                           for(int a=numInterSet; a>k; a--)
                           {
                              interSetStart[a]=interSetStart[a]+1;
                           }

                           for(int a=k+1; a<numInterSet; a++)
                           {
                              int *temp=validIntersections[interSetStart[a]-1];
                              // printf("interSetStart[a]-1 = %d",interSetStart[a]-1);
                              validIntersections[interSetStart[a]-1] = validIntersections[numIntersection];
                              validIntersections[numIntersection] = temp;
                           }
                           break;
                        }
                        
                        //Add new set at the end
                        if(k==numInterSet-1 && intersectionLen < interSetLen[k])
                        {

                           //printf("Case 4 new set at the end\n");
                           interSetLen[numInterSet] = intersectionLen;
                           interSetStart[numInterSet] = numIntersection;
                           numInterSet++;
                           break;
                        }
                     }
                  }
                  numIntersection++;
                  // printf("Updated len and start:\n");
                  // print_array_int(interSetLen, numInterSet);
                  // print_array_int(interSetStart, numInterSet);
                  // print_array_int_int(validIntersections, numIntersection, interSetStart, interSetLen, numInterSet);
                  // printf("Total length %d\n", numIntersection);
                  // printf("breakpoint\n");
                  ;
               } 
               else
               {
                  printf("free\n");
                  free(validIntersections[numIntersection]);
               }
            }  
         }
         print_array_int(interSetLen, numInterSet);
         print_array_int(interSetStart, numInterSet);
         //print_array_int_int(validIntersections, numIntersection, interSetStart, interSetLen, numInterSet);
         
         printf("Total length %d max length %d\n", numIntersection, inst->num_rows*(inst->num_rows-1)/2);

         FILE *f;
         f = fopen("Intersections.txt", "w");

         fprint_array_int_int(f, validIntersections, numIntersection, interSetStart, interSetLen, numInterSet);
         fclose(f);               
         inst->intersections = validIntersections;
         inst->interSetLen = interSetLen;
         inst->interSetStart = interSetStart;
         inst->numInterSet = numInterSet;
	      inst->numIntersections = numIntersection;


         printf("Installing generic callback\n");
         CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_BRANCHING, genericcallbackfunc, inst);
         int ncores = 1;
         //CPXgetnumcores(env, &ncores);
         CPXsetintparam(env, CPX_PARAM_THREADS, ncores); // it was reset after callback
   
      }

      // CPXwriteprob(env, presolve, "modelpresolved.lp", NULL);
      if (CPXmipopt(env, lp)) // solve the problem
         printf(" Problems on CPXmipopt");

      CPXwriteprob(env, lp, "model.lp", NULL);   
   }
   


   // free and close CPLEX model
   CPXfreeprob(env, &lp);
   CPXcloseCPLEX(&env);

}


/** 
 * Callback function only useful save the preprocessed redced problem as a new one.
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
 * @param useraction_p 
 * 
 * 
 */
int preprocessinglegacycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                         const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
   {
      instance *inst = (instance *)cbhandle; // casting of cbhandle
      // printf("branching type = %c\n", (char)type);
      // printf("cplex wants to create %d new branches\n", nodecnt);
      // printf("in total changes %d\n", bdcnt);
      // printf("cplex chose the variable: %d, to be bound\n", indices[0]);
      // printf("in the first node it is set to %f \n", bd[0]);
      // printf("in the second node it is set to %f \n", bd[1]);
      int seqnum_p = 0;

      int indice1 = indices[0];
      int indice2 = indices[0];
      char lu1 = lu[0];
      char lu2 = lu[1];
      double bd1 = bd[0];
      double bd2 = bd[1];

      CPXLPptr nodelp;
      CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);
      //int  CPXgetrows( CPXCENVptr env, CPXCLPptr lp, int * nzcnt_p, int * rmatbeg, int * rmatind, double * rmatval, int rmatspace, int * surplus_p, int begin, int end )
      int error;
      //CPXLPptr clone = CPXcloneprob( env, nodelp, &error);
      CPXLPptr clone = nodelp;
      int num_cols = CPXgetnumcols(env, clone);
      int varInd[num_cols];
      char ctype[num_cols];
      for(int i=0; i<num_cols; i++)
      {
         varInd[i]=i;
         ctype[i]='B';
      }
      int currthread = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &currthread);  
      int solvednodes = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &solvednodes);  
      //int iterations = -1; CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MIP_ITERATIONS, &iterations);  
      
      CPXchgctype (env, clone, num_cols, varInd, ctype);
      int num_rows = CPXgetnumrows(env, clone);
      if(inst->num_cols != num_cols || inst->num_rows != num_rows)
      {
         printf("Thread %d. Solved nodes %d. The presolve performed a change in the problem now there are %d columns and %d rows\n", currthread, solvednodes, num_cols, num_rows);
         inst->num_cols=num_cols;
         inst->num_rows=num_rows;
      }
      else
      {
         printf("Thread %d. Solved nodes %d. No changes to the model\n", currthread, solvednodes);
      }
      // printf("the presolved problem has %d columns\n", num_cols);
      // printf("the presolved problem has %d rows\n", num_rows);
      char newname[100];
      sprintf(newname, "../Dataset/ReducedLP%.*s_reduced.lp", (int)(strlen(inst->input_file)-3-13), inst->input_file+13);
      printf("%s\n", inst->input_file);
      printf("%s\n", newname);
      
      CPXwriteprob(env, clone, newname, NULL);
      printf("Problem printed\n");
      
      return 0;   
   }
