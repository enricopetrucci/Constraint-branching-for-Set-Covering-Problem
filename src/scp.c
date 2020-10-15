/**
 * This file contains functions needed for computing the solutions to the
 * scp using the CPLEX library.
 * 
 * @author Petrucci Enrico
*/


#include "scp.h"

/**
 * It specifies the behavior for the generic callback in the branching context, that is once the LP relaxation 
 * is solved and CPLEX decided that it is time to branch.
 *  
 * @param context CPLEX context in which the callback has been called
 * @param contextid CPLEX context id
 * @param userhandle pointer to our own struct inst that we gave CPLEX upon linking out callback function
 */
int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle) 
{
   if ( contextid == CPX_CALLBACKCONTEXT_BRANCHING)
   {
      // cast cbhandle to our structure of type instance
      instance *inst = (instance *)cbhandle; // casting of cbhandle
         
      if(inst->branching==0)
      {
         return 0;
      }
      else
      {
         // cast cbhandle to our structure of type instance
         instance *inst = (instance *)cbhandle; // casting of cbhandle
         
         /* Check the status of the continuous relaxation. */
         int statind;
         // with 0 it forces the optimality of the lp-relaxation
         CPXcallbackgetrelaxationstatus(context, &statind, 0);
         int threadNum;
         CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &threadNum);

         if (statind == CPX_STAT_OPTIMAL || statind == CPX_STAT_OPTIMAL_INFEAS)
         {
            // Get LP relaxation solution and obj
            double *x = (double *) calloc(inst->num_cols, sizeof(double));
            double obj;
            if(CPXcallbackgetrelaxationpoint(context, x, 0, inst->num_cols-1, &obj)!=0)
               print_error("Problems in CPXcallbackgetrelaxationpoint\n");
            //printf("LP relaxation solved optimally it has the objective %f\n", obj);
            //print_array(x, inst->num_cols);
   
            //perform constraint branching
            if(inst->branching==1)
            {
               int foundConstraint = 0;
               int i = 0;
               int j = 0;

               switch(inst->constraintBranchVer)
               {
                  case 0: //always use constraint branching on the longest violated intersection.
                     findBranchingConstraint0(&i, &j, x, inst);
                     break;
                  case 1:
                     findBranchingConstraint1(&i, &j, x, inst, 0.005);
                     break;
                  case 2:
                  {
                     double *ub = (double *) calloc(inst->num_cols, sizeof(double));
                     double *lb = (double *) calloc(inst->num_cols, sizeof(double));

                     CPXcallbackgetlocalub(context, ub, 0, inst->num_cols-1);
                     CPXcallbackgetgloballb (context, lb, 0, inst->num_cols-1);

                     int *fixed = (int *) calloc(inst->num_cols, sizeof(int));
                     int count = 0;
                     for(int i=0; i<inst->num_cols; i++)
                     {
                        if(ub[i]==lb[i])
                           {
                              fixed[i]=1;
                              count++;
                           }
                        else
                           fixed[i]=0;
                     }
                     //printf("locally %d variables are fixed out of %d\n", count, inst->num_cols);
                     findBranchingConstraint2(&i, &j, x, inst, fixed);
                     break;
                  }
                  case 3:
                     findBranchingConstraint3(&i, &j, x, inst);
                     break;
                  case 4:
                  {
                     double *ub = (double *) calloc(inst->num_cols, sizeof(double));
                     double *lb = (double *) calloc(inst->num_cols, sizeof(double));

                     CPXcallbackgetlocalub(context, ub, 0, inst->num_cols-1);
                     CPXcallbackgetgloballb (context, lb, 0, inst->num_cols-1);

                     int *fixed = (int *) calloc(inst->num_cols, sizeof(int));
                     int count = 0;
                     for(int i=0; i<inst->num_cols; i++)
                     {
                        if(ub[i]==lb[i])
                           {
                              fixed[i]=1;
                              count++;
                           }
                        else
                           fixed[i]=0;
                     }
                     //printf("locally %d variables are fixed out of %d\n", count, inst->num_cols);
                     findBranchingConstraint4(&i, &j, x, inst, fixed);
                     break;
                  }
                  case 5:
                  {
                     double *ub = (double *) calloc(inst->num_cols, sizeof(double));
                     double *lb = (double *) calloc(inst->num_cols, sizeof(double));

                     CPXcallbackgetlocalub(context, ub, 0, inst->num_cols-1);
                     CPXcallbackgetgloballb (context, lb, 0, inst->num_cols-1);

                     int *fixed = (int *) calloc(inst->num_cols, sizeof(int));
                     int count = 0;
                     for(int i=0; i<inst->num_cols; i++)
                     {
                        if(ub[i]==lb[i])
                           {
                              fixed[i]=1;
                              count++;
                           }
                        else
                           fixed[i]=0;
                     }
                     //printf("locally %d variables are fixed out of %d\n", count, inst->num_cols);
                     findBranchingConstraint5(&i, &j, x, inst, fixed);
                     break;       
                  }           
               }
               
               if (i!= -1)
               {
                  inst->constraintBranching[threadNum]++;
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
                  if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs1, "E", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child1)!=0)
                  {
                     print_error("First branch was not succesful\n");
                  }
                  //if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs2, "E", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child2)==0)
                  if(CPXcallbackmakebranch(context, 0, NULL, NULL, NULL, 1, inst->interSetLen[i], &rhs2, "G", &izero, inst->intersections[inst->interSetStart[i]+j], values, obj, &child2)!=0)
                  {
                     print_error("Second branch was not succesful\n");
                  }
                  free(values);
               }
               else
               {
                  inst->defaultBranching[threadNum]++;
               }
            }
            
            else if(inst->branching==-1)
            // perform most infeasible branching
            {
               double maxfrac = 0.0;
               int maxvar = -1;

               for(int i=0; i<inst->num_cols; i++)
               {
                  double intval = round (x[i]);
                  double frac = fabs (intval - x[i]);

                  if ( frac > maxfrac ) 
                  {
                     maxfrac = frac;
                     maxvar = i;
                  }
               }
               // printf("fixing variable %d that has value %f\n", maxvar, x[maxvar]);
               int child1, child2;
               double const up = ceil (x[maxvar]);
               double const down = floor (x[maxvar]);
               CPXcallbackmakebranch(context, 1, &maxvar, "L", &up, 0, 0, NULL, NULL, NULL, NULL, NULL, obj, &child1);
               CPXcallbackmakebranch(context, 1, &maxvar, "U", &down, 0, 0, NULL, NULL, NULL, NULL, NULL, obj, &child2);
            }
            free(x);
         }
         else 
         {
            printf("LP relaxation not solved optimally\n");
         }
      }
      //printf("callback END\n");
      return 0;
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
   char str[80];

   if(inst->branching==0)
   {
      sprintf(str, "logfile_defaultBranching.txt");
   }
   else
   {
      sprintf(str, "logfile_constraintBranching_version%d.txt", inst->constraintBranchVer);
   }
   
   CPXsetlogfilename(env, str, "w");

   // import problem from file
   CPXreadcopyprob(env, lp, inst->input_file, "LP");
   
   inst->num_cols = CPXgetnumcols(env, lp);
   inst->num_rows = CPXgetnumrows(env, lp);
   
   if(inst->extractPreprocessing==1)
   {
      performPreprocessing(env, lp, inst);   
   }
   else
   {
      printf("Installing generic callback\n");
      CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_BRANCHING, genericcallbackfunc, inst);
            
      if (inst->threads==0)
         CPXgetnumcores(env, &inst->threads);

      CPXsetintparam(env, CPX_PARAM_THREADS, inst->threads);

      inst->constraintBranching=(int *) calloc(inst->threads, sizeof(int));
      inst->defaultBranching=(int *) calloc(inst->threads, sizeof(int));


      if(inst->branching == 1)
      {

         // Get the rows of the problem
         int nnz;
         int *izero = (int *) calloc(inst->num_rows, sizeof(int));
         int *indexes = (int *) calloc(inst->num_rows*inst->num_cols, sizeof(int));// this one is big!
         double *values = (double *) calloc(inst->num_rows*inst->num_cols, sizeof(double));// this one is big too!
         int space;
         int surplus_p;

         CPXgetrows(env, lp, &nnz, izero, indexes, values, inst->num_rows*inst->num_cols, &surplus_p, 0, inst->num_rows-1);

         indexes = (int *) realloc(indexes, nnz*sizeof(int));// realloc should free the unused spaced of the array
         values = (double *) realloc(values, nnz*sizeof(double));// realloc should free the unused spaced of the array
         double start = second();
         populateIntesectionsOf2(izero, indexes, nnz, inst);
         double end = second();
         printf("Found %d constraint intersections in %f\n", inst->numIntersections, end-start);
         
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
         
         print_array_int(inst->variableFreq, inst->num_cols);
         // int ncores = 1;
         // //CPXgetnumcores(env, &ncores);
         // CPXsetintparam(env, CPX_PARAM_THREADS, ncores); // it was reset after callback
         
         free(izero);
         free(indexes);
         free(values);
      
      }
   
      if (CPXmipopt(env, lp)) // solve the problem
         print_error(" Problems on CPXmipopt");

      //CPXwriteprob(env, lp, "model.lp", NULL);   
      inst->solution = (double *) calloc(inst->num_cols, sizeof(double));
      double obj;
         
      if(CPXgetx(env, lp, inst->solution, 0, inst->num_cols-1)==0 && CPXgetbestobjval(env, lp, &obj)==0)
         printf("Solution found with objval = %f\n", obj);
      else
         print_error("Problems in getting the solution");
      
      if(inst->branching == 1)
      {
         int constraintCount=0;
         int defaultCount=0;

         for(int i=0; i<inst->threads; i++)
         {
            constraintCount += (inst->constraintBranching[i]);
            defaultCount += (inst->defaultBranching[i]);
         }
         printf("%d constraint branching\n%d default branching\n", constraintCount, defaultCount);
         for(int i=0; i<inst->numIntersections; i++)
         {
            free(inst->intersections[i]);      
         }
         free(inst->intersections);
         free(inst->interSetLen);
         free(inst->interSetStart);
      }
      free(inst->solution);
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
      
      return 1;   
   }



/** 
 * It performes and save the preprocessing of a problem.
 * All the cuts generated by CPLEX are disabled.
 * 
 *  
 * @param env environment
 * @param lp problem
 * @param inst instance of our problem
 * 
 */
void performPreprocessing(CPXENVptr env, CPXLPptr lp, instance* inst)
{
   //install legacy callback for extracting the preprocessed problem
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

   printf("The original models has %d columns\n", inst->num_cols);
   printf("The original models has %d rows\n", inst->num_rows);

   if (CPXmipopt(env, lp)) // solve the problem
      print_error(" Problems on CPXmipopt");
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




void findBranchingConstraint0(int* n, int* m, double* x, instance* inst)
{
   int foundConstraint = 0;
   int i = 0;
   int j = 0;
   
   // cycle on all the sets
   for(; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      j=0;
      // cycle on all the intersection having the same lenght
      for (; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         //if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i] || sum>1+1e-05*inst->interSetLen[i]) // set partitionoing case
         if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i]) // set covering
         {
            //printf("Found constraint not yet used: i=%d, j=%d\n", i, j);
            foundConstraint=1;
            break;
         }
      }
      if(foundConstraint==1)
      {
         break;
      }
   }
   if(foundConstraint == 1)
   {
      *n=i;
      *m=j;
   }
   
}



/**
 * Branch on the constraint with the highest number of variables that are in the interval 0.5+-delta
 * 
 * 
 */

void findBranchingConstraint1(int* n, int* m, double* x, instance* inst, double delta)
{
   int foundConstraint = 0;
   
   int best_i=-1;
   int best_j=-1;
   int MaxHighInfeas = 0;


   // cycle on all the sets
   for(int i = 0; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      // cycle on all the intersection having the same lenght
      for (int j = 0; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         int countHighInfeas = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
            if(x[inst->intersections[inst->interSetStart[i]+j][k]]>(0.5-delta) && x[inst->intersections[inst->interSetStart[i]+j][k]]<0.5+delta)
            {
               countHighInfeas++;
            }
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i] && countHighInfeas>MaxHighInfeas) // the constraint is violated
         {
            //printf("Found constraint not yet used: i=%d, j=%d, countHighInfeas = %d\n", i, j, countHighInfeas);
            best_i=i;
            best_j=j;
            MaxHighInfeas = countHighInfeas;
            foundConstraint++;
         }
      }
   }
   if(foundConstraint > 0)
   {
      //printf("constraint updated %d times\n", foundConstraint);
      *n = best_i;
      *m = best_j;
   }
   else
   {
      //printf("performing standard branching\n");
      *n=-1;
      *m=-1;
   }
   //printf("breakpoint\n");
   
}





/**
 * Branch on the constraint with the highest number of unfixed variables
 * 
 * 
 */

void findBranchingConstraint2(int* n, int* m, double* x, instance* inst, int* fixed)
{
   double delta=0;
   int foundConstraint = 0;
   
   int best_i=-1;
   int best_j=-1;
   int MaxFree = 0;


   // cycle on all the sets
   for(int i = 0; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      // cycle on all the intersection having the same lenght
      for (int j = 0; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         int countFixed = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
            countFixed += fixed[inst->intersections[inst->interSetStart[i]+j][k]];
            
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i] && (inst->interSetLen[i]-countFixed) > MaxFree) // the constraint is violated
         {
            //printf("Found constraint not yet used: i=%d, j=%d, countHighInfeas = %d\n", i, j, countHighInfeas);
            best_i=i;
            best_j=j;
            MaxFree = inst->interSetLen[i]-countFixed;
            foundConstraint++;
         }
      }
   }
   if(foundConstraint > 0)
   {
      if( MaxFree < inst->lowestNumVariables)
      {
         inst->lowestNumVariables = MaxFree;
         printf("branching on a constraint with %d unfixed variables\n", MaxFree);
         
      }
      //printf("constraint updated %d times: branching constraint has %d unfixed variables out of the total %d\n", foundConstraint, MaxFree, inst->interSetLen[best_i]);
      *n = best_i;
      *m = best_j;
   }
   else
   {
      //printf("performing standard branching\n");
      *n=-1;
      *m=-1;
   }
   //printf("breakpoint\n");
   
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
 * Compute all possible intersections between the 
 * constraint and stored them sorted by length inside inst->intersections
 *  
 * @param izero starting index for each constraint
 * @param indexes indexes of the variables in the constraints
 * @param nnz number of non zeros inside index
 * @param inst instance of our problem
 * 
 */   
void populateIntesectionsOf4(instance* inst)
{
   int maxlen = INT_MAX;//inst->numIntersections*(inst->numIntersections-1)/2;
   printf("maxlen = %d\n", maxlen);
   int** intersections = (int **)calloc(maxlen, sizeof(int*)); 
   int* interSetLen = (int *)calloc(inst->num_cols, sizeof(int));
   int* interSetStart = (int *)calloc(inst->num_cols, sizeof(int));
   int numInterSet = 0;
   int numIntersections = 0;
      
   int index_set_first = 0;
   int index_set_second = 0;
   
   // cycle the constraints, choosing the first one
   for(int i = 0; i < (inst -> numIntersections-1); i++)
   {
      if(i==inst->interSetStart[index_set_first+1])
         index_set_first++;
      //printf("Considering i = %d\n", i);
      //Select the second constraint
      for(int n = i+1; n < inst -> numIntersections; n++)
      {
         //printf("Considering n = %d\n", n);
         if(n==inst->interSetStart[index_set_second+1])
         index_set_second++;
         
         int j = 0;
         int m = 0;
         int intersectionLen = 0;

         //printf("Compared with n = %d, that starts in %d it has length %d\n", n, izero[n], length);
         
         // alloc space for current intersection
         intersections[numIntersections] = (int *) calloc(inst->num_cols, sizeof(int));
         
         // cycle on all the variables inside both constraints
         while(j < inst->interSetLen[index_set_first] && m < inst->interSetLen[index_set_second])// problem for the last constraint's length
         {  
            //printf("Considering j = %d, m = %d\n", j, m);
            // case 1: same variable. The variable is added to the computed information
            if(inst->intersections[i][j]==inst->intersections[n][m])
            {
               // printf("found intersection for variable %d in positions %d and %d\n", indexes[izero[i]+j], j, m);
               intersections[numIntersections][intersectionLen]=inst->intersections[n][m];
               intersectionLen++;
               m++;
               j++;
            }
            // case 2 and 3: mismatch, increase the index of corresponding to the lowest variable
            if(inst->intersections[i][j]<inst->intersections[n][m])
            {
               j++;
            }
            
            if(inst->intersections[i][j]>inst->intersections[n][m])
            {
               m++;
            }
         }
         // if the intersecion is non empty sort the validIntersection array and update the other arrays accordingly
         if(intersectionLen>1)
         {
            // realloc should free the unused spaced of the array
            intersections[numIntersections] = (int *) realloc(intersections[numIntersections], intersectionLen*sizeof(int));
            // printf("Intersection between constraints %d and %d has %d variables\n", i, n, intersectionLen);
            // case 1 validIntersection is empty
            if(numIntersections==0)
            {
               // printf("Case 1 empty list\n");
               interSetLen[0]=intersectionLen;
               interSetStart[0]=0;
               numInterSet++;
            }
            else
            {
               // cycle on all sets
               for(int k=0; k < numInterSet; k++)
               {
                  // Case 2 Add new set at the beginning or add new set in between two sets
                  if((k==0 && intersectionLen > interSetLen[k]) || (k>0 && intersectionLen<interSetLen[k-1] && intersectionLen>interSetLen[k]))
                  {
                     // printf("Case 2 new set in between\n");

                     // Update information on inst->interSetLen and inst->interSetStart
                     for(int a=numInterSet; a>=k; a--)
                     {
                        interSetLen[a+1]=interSetLen[a];
                        interSetStart[a+1]=interSetStart[a]+1;
                     }

                     interSetLen[k]=intersectionLen;
                     numInterSet++;
            
                     // Updae infromation in validIntersection on the right with respect to k
                     for(int a=k; a<numInterSet-1; a++)
                     {
                        // printf("from last position to %d\n", inst->interSetStart[a+1]-1);
                        int *temp=intersections[interSetStart[a+1]-1];
                        intersections[interSetStart[a+1]-1] = intersections[numIntersections];
                        intersections[numIntersections] = temp;
                     }
                     break;
                  }
                  
                  // Case 3 no new set needed, only add an intersection and update everithing accordingly
                  if(intersectionLen == interSetLen[k])
                  {

                     // printf("Case 3 same set\n");
                     for(int a=numInterSet; a>k; a--)
                     {
                        interSetStart[a]=interSetStart[a]+1;
                     }

                     for(int a=k+1; a<numInterSet; a++)
                     {
                        int *temp=intersections[interSetStart[a]-1];
                        // printf("inst->interSetStart[a]-1 = %d",inst->interSetStart[a]-1);
                        intersections[interSetStart[a]-1] = intersections[numIntersections];
                        intersections[numIntersections] = temp;
                     }
                     break;
                  }
                  
                  // Case 4 Add new set at the end
                  if(k==numInterSet-1 && intersectionLen < interSetLen[k])
                  {

                     //printf("Case 4 new set at the end\n");
                     interSetLen[numInterSet] = intersectionLen;
                     interSetStart[numInterSet] = numIntersections;
                     numInterSet++;
                     break;
                  }
               }
            }
            numIntersections++;
            //printf("%d\n",numIntersections);
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
            free(intersections[numIntersections]);
         }
      }  
   }
   printf("Length for each set");
   print_array_int(interSetLen, numInterSet);
   printf("Starting index for each set for each set");
   print_array_int(interSetStart, numInterSet);
   //print_array_int_int(inst->intersections, inst->numIntersections, inst->interSetStart, inst->interSetLen, inst->numInterSet);
   
   printf("Total length %d max length %d\n", numIntersections, inst->numIntersections*(inst->numIntersections-1)/2);

   FILE *f;
   f = fopen("Intersections4.txt", "w");

   fprint_array_int_int(f, intersections, numIntersections, interSetStart, interSetLen, numInterSet);
   fclose(f);
   inst->shortestConstraint = interSetLen[0];
   inst->lowestNumVariables = interSetLen[0];


   inst->intersections = intersections; 
   inst->interSetLen = interSetLen;
   inst->interSetStart = interSetStart;
   inst->numInterSet = numInterSet;
   inst->numIntersections = numIntersections;
}




/**
 * Branch on the constraint with the highest number of variables that are in the interval 0.5+-delta
 * 
 * 
 */
void findBranchingConstraint3(int* n, int* m, double* x, instance* inst)
{
   int foundConstraint = 0;
   
   int best_i=-1;
   int best_j=-1;
   double MinIntersectionValue = 1.0;


   // cycle on all the sets
   for(int i = 0; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      // cycle on all the intersection having the same lenght
      for (int j = 0; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         // printf("Comparing lowest sum %f with current sum %f\n", MinIntersectionValue, sum);
         if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i] && sum<MinIntersectionValue) // the constraint is violated
         {
            //printf("Found constraint not yet used: i=%d, j=%d, countHighInfeas = %d\n", i, j, countHighInfeas);
            best_i=i;
            best_j=j;
            MinIntersectionValue = sum;
            foundConstraint++;
         }
      }
   }
   if(foundConstraint > 0)
   {
      if(foundConstraint>1)
         // printf("constraint updated %d times\n", foundConstraint);
      *n = best_i;
      *m = best_j;
   }
   else
   {
      //printf("performing standard branching\n");
      *n=-1;
      *m=-1;
   }
   // printf("breakpoint\n");
}



/**
 * Branch on the constraint with the highest number of unfixed variables
 * 
 * 
 */

void findBranchingConstraint4(int* n, int* m, double* x, instance* inst, int* fixed)
{
   double delta=0;
   int foundConstraint = 0;
   
   int best_i=-1;
   int best_j=-1;
   int MaxFree = 0;
   double MinIntersectionValue = 1.0;



   // cycle on all the sets
   for(int i = 0; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      // cycle on all the intersection having the same lenght
      for (int j = 0; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         int countFixed = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
            countFixed += fixed[inst->intersections[inst->interSetStart[i]+j][k]];
            
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         if(sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i] && ((inst->interSetLen[i]-countFixed > MaxFree) || (inst->interSetLen[i]-countFixed == MaxFree && sum<MinIntersectionValue))) // the constraint is violated
         {
            //printf("Found constraint not yet used: i=%d, j=%d, countHighInfeas = %d\n", i, j, countHighInfeas);
            best_i=i;
            best_j=j;
            MaxFree = inst->interSetLen[i]-countFixed;
            MinIntersectionValue = sum;
       
            foundConstraint++;
         }
      }
   }
   if(foundConstraint > 0)
   {
      if( MaxFree < inst->lowestNumVariables)
      {
         inst->lowestNumVariables = MaxFree;
         printf("branching on a constraint with %d unfixed variables\n", MaxFree);
         
      }
      //printf("constraint updated %d times: branching constraint has %d unfixed variables out of the total %d\n", foundConstraint, MaxFree, inst->interSetLen[best_i]);
      *n = best_i;
      *m = best_j;
   }
   else
   {
      //printf("performing standard branching\n");
      *n=-1;
      *m=-1;
   }
   //printf("breakpoint\n");  
}



/**
 * Branch on the constraint with the highest number of unfixed variables
 * 
 * 
 */

void findBranchingConstraint5(int* n, int* m, double* x, instance* inst, int* fixed)
{
   double delta=0;
   int foundConstraint = 0;
   
   int best_i=-1;
   int best_j=-1;
   int MaxFree = 0;
   double MinIntersectionValue = 1.0;

   int mostFrequentUnfixedVar = 0;

   for(int i=0; i<inst->num_cols; i++)
   {
      if(fixed[inst->variableFreq[i]]==0)
      {
         mostFrequentUnfixedVar = i;
         break;
      }
   }


   // cycle on all the sets
   for(int i = 0; i < inst->numInterSet; i++)
   {
      // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
      // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
      // cycle on all the intersection having the same lenght
      for (int j = 0; j<( (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]); j++)
      {
         // printf("j=%d\n", j);
         double sum = 0;
         int containsFreqVar = 0;
         // cycle all the variables in the intersection and get the sum
         for(int k = 0; k < inst->interSetLen[i]; k++)
         {
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=x[inst->intersections[inst->interSetStart[i]+j][k]];
            if(inst->intersections[inst->interSetStart[i]+j][k]==mostFrequentUnfixedVar)
            {
               containsFreqVar = 1;
            }
         }
         // printf("Constraint: i=%d, j=%d has sum in the solution = %.20f \n", i, j, sum);   
         
         // if the sum of the variables in the intersection is strictly (considering tolerance) in between 0 and 1 the constraint can be added.
         if(containsFreqVar==1 && sum>1e-05*inst->interSetLen[i] && sum<1-1e-05*inst->interSetLen[i]) // the constraint is violated
         {
            //printf("Found constraint not yet used: i=%d, j=%d, countHighInfeas = %d\n", i, j, countHighInfeas);
            best_i=i;
            best_j=j;
            foundConstraint++;
            break;
         }
      }
   }
   if(foundConstraint > 0)
   {
      if( MaxFree < inst->lowestNumVariables)
      {
         inst->lowestNumVariables = MaxFree;
         printf("branching on a constraint with %d unfixed variables\n", MaxFree);
         
      }
      //printf("constraint updated %d times: branching constraint has %d unfixed variables out of the total %d\n", foundConstraint, MaxFree, inst->interSetLen[best_i]);
      *n = best_i;
      *m = best_j;
   }
   else
   {
      //printf("performing standard branching\n");
      *n=-1;
      *m=-1;
   }
   //printf("breakpoint\n");  
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
void computeVariableFrequency(int* izero, int* indexes, int nnz, instance* inst)
{
   int* frequencies = (int *)calloc(inst->num_cols, sizeof(int));
   int* variables = (int *)calloc(inst->num_cols, sizeof(int));
   
   // cycle the constraints, choosing the first one
   for(int i = 0; i < nnz; i++)
   {
      frequencies[indexes[i]]++;
   }
   
   print_array_int(frequencies, inst->num_cols);
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
}