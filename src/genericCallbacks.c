/**
 * This file contains functions that are used as CPLEX generic 
 * callbacks and other related functions.
 *
 * @author Petrucci Enrico
*/


#include "scp.h"



void solveUsingGenericCallback(CPXENVptr env, CPXLPptr lp, instance* inst)
{
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
      int *indexes = (int *) calloc(inst->num_rows*inst->num_cols, sizeof(int));
      double *values = (double *) calloc(inst->num_rows*inst->num_cols, sizeof(double));
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
      end = second();
      printf("Computed variable frequencies in %f\n", end-start);

      free(izero);
      free(indexes);
      free(values);
   }

   if(CPXmipopt(env, lp)) // solve the problem
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
         for(int i=0; i<inst->numIntersections; i++)
         {
            free(inst->intersections[i]);
         }
         free(inst->intersections);
         free(inst->interSetLen);
         free(inst->interSetStart);
         free(inst->variableFreq);
   }
   int constraintCount=0;
   int defaultCount=0;

   for(int i=0; i<inst->threads; i++)
   {
         constraintCount += (inst->constraintBranching[i]);
         defaultCount += (inst->defaultBranching[i]);
   }
   printf("%d constraint branching\n%d default branching\n", constraintCount, defaultCount);

   free(inst->solution);
   free(inst->constraintBranching);
   free(inst->defaultBranching);
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
 * Branch on the constraint with highest number of unfixed variables and lowest sum value
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
 * Branch on the constraint with the highest number of variables that contains
 * the most frequent unfixed variable
 *
 */

void findBranchingConstraint5(int* n, int* m, double* x, instance* inst, int* fixed)
{
   int foundConstraint = 0;

   int best_i=-1;
   int best_j=-1;

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
 * It specifies the behavior for the generic callback in the branching context, that is once the LP relaxation
 * is solved and CPLEX decided that it is time to branch.
 *
 * @param context CPLEX context in which the callback has been called
 * @param contextid CPLEX context id
 * @param userhandle pointer to our own struct inst that we gave CPLEX upon linking out callback function
 */
int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle)
{
   if (contextid == CPX_CALLBACKCONTEXT_BRANCHING)
   {
      // cast cbhandle to our structure of type instance
      instance *inst = (instance *)cbhandle; // casting of cbhandle

      if(inst->branching==0)
      {
         int threadNum;
         CPXcallbackgetinfoint(context, CPXCALLBACKINFO_THREADID, &threadNum);
         inst->defaultBranching[threadNum]++;
         return 0;
      }
      else
      {
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

                     free(ub);
                     free(lb);
                     free(fixed);

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

                     free(ub);
                     free(lb);
                     free(fixed);
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

                     free(ub);
                     free(lb);
                     free(fixed);

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
