/**
 * This file contains functions needed for computing the solutions to the
 * scp using the CPLEX library.
 * 
 * @author Petrucci Enrico
*/


#include "scp.h"

int genericcallbackfunc(CPXCALLBACKCONTEXTptr context, CPXLONG contextid, void *cbhandle) 
{
    if ( contextid == CPX_CALLBACKCONTEXT_BRANCHING )
    {
       printf("inside branching generic callback\n");
      //   /* Check the status of the continuous relaxation. */
      //   int statind;

      //   // with 0 it forces the optimality of the lp-relaxation
      //   CPXcallbackgetrelaxationstatus(context, &statind, 0);
      //   if ( statind == CPX_STAT_OPTIMAL || statind == CPX_STAT_OPTIMAL_INFEAS   )
      //   {
      //    /* Continuous relaxation is solved to optimality. Get the current
      //     * x-vector and the objective value and create custom branches.
      //     */
      //    double *x = ;
      //    double obj;
      //    CPXcallbackgetrelaxationpoint(context, x, 0, nvars - 1, &obj);
      //    /*
      //    probably never needed
      //    if ( ... ) {
      //       /* We decided to cut off the current node. */
      //       CPXcallbackrelaxationprunenode(context);
      //    /*
      //    }
      //    else {
      //       /* We want to branch. */
      //       CPXCNT child1, child2;
      //       CPXcallbackmakebranch(context, ..., obj, &child1);
      //       CPXcallbackmakebranch(context, ..., obj, &child2);
      //    }
      // }
      // else {
      //    /* The continuous relaxation was not solved to optimality.
      //     * There may be a problem and some special handling may be required.
      //     */
      //    ...
      // }
   }
}


/** 
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
int branchlegacycallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
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
      CPXLPptr clone = CPXcloneprob( env, nodelp, &error);
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
         printf("Thread %d. Branching %d. Solved nodes %d. The presolve performed a change in the problem now there are %d columns and %d rows\n", currthread, inst->branchcalls[currthread]++, solvednodes, num_cols, num_rows);
         inst->num_cols=num_cols;
         inst->num_rows=num_rows;
      }
      else
      {
         printf("Thread %d. Branching %d. Solved nodes %d. No changes to the model\n", currthread, inst->branchcalls[currthread]++, solvednodes); 
      }
      // printf("the presolved problem has %d columns\n", num_cols);
      // printf("the presolved problem has %d rows\n", num_rows);
      // CPXwriteprob(env, clone, "nodelp.lp", NULL);
      // printf("Problem printed\n");
      
      double rhs = 0;
      char sense = 'L';
      int rmatbeg = 0;
      int rmatind = 0;
      double rmatval = 0;
      //CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, 1, &indice1, &lu1, &bd1, 0, 0, &rhs, &sense, &rmatbeg, &rmatind, &rmatval, nodeest[0], NULL, &seqnum_p);
      //CPXbranchcallbackbranchgeneral(env, cbdata, wherefrom, 1, &indice2, &lu2, &bd2, 0, 0, &rhs, &sense, &rmatbeg, &rmatind, &rmatval, nodeest[1], NULL, &seqnum_p);
      
      CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, &indice1, &lu1, &bd1, nodeest[0], NULL, &seqnum_p);
      CPXbranchcallbackbranchbds(env, cbdata, wherefrom, 1, &indice2, &lu2, &bd2, nodeest[1], NULL, &seqnum_p);
      //*useraction_p = 0;
   }


/** 
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
         printf("Thread %d. Branching %d. Solved nodes %d. The presolve performed a change in the problem now there are %d columns and %d rows\n", currthread, inst->branchcalls[currthread]++, solvednodes, num_cols, num_rows);
         inst->num_cols=num_cols;
         inst->num_rows=num_rows;
      }
      else
      {
         printf("Thread %d. Branching %d. Solved nodes %d. No changes to the model\n", currthread, inst->branchcalls[currthread]++, solvednodes); 
      }
      // printf("the presolved problem has %d columns\n", num_cols);
      // printf("the presolved problem has %d rows\n", num_rows);
      char newname[100];
      sprintf(newname, "../Dataset/ReducedLP%.*s_reduced.lp", (strlen(inst->input_file)-3-13), inst->input_file+13);
      printf("%s\n", inst->input_file);
      printf("%s\n", newname);
      
      CPXwriteprob(env, clone, newname, NULL);
      printf("Problem printed\n");
      
      return 0;   
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
   inst->branchcalls[0]=0;
   inst->branchcalls[1]=0;

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
         printf("Installing legacy callback\n");
         //CPXsetintparam(env, CPX_PARAM_MIPCBREDLP, CPX_OFF); // let MIP callbacks work on the original model
         CPXsetintparam(env, CPX_PARAM_REPEATPRESOLVE, 1); // let MIP callbacks work on the original model
         
         CPXsetbranchcallbackfunc(env, branchlegacycallback, inst );
         int ncores = 2;
         // CPXgetnumcores(env, &ncores);
         CPXsetintparam(env, CPX_PARAM_THREADS, ncores); // it was reset after callback
      }
      else if(inst->callback==2)
      {
         printf("Installing generic callback\n");
         // probably never needed
      
         CPXcallbacksetfunc(env, lp, CPX_CALLBACKCONTEXT_BRANCHING, genericcallbackfunc, inst);
      }

      // CPXLPptr lpPresolve = CPXcreateprob(env, &error, "scppresolve");
      // CPXreadcopyprob(env, lpPresolve, inst->input_file, "LP");
      // CPXpresolve(env, lpPresolve, CPX_ALG_NONE);
      // CPXLPptr presolve = CPXcreateprob(env, &error, "presolve");
      // CPXgetredlp( env, lpPresolve,  &presolve);

      // CPXwriteprob(env, presolve, "modelpresolved.lp", NULL);
      if (CPXmipopt(env, lp)) // solve the problem
         printf(" Problems on CPXmipopt");

      CPXwriteprob(env, lp, "model.lp", NULL);   
   }
   


   // free and close CPLEX model
   CPXfreeprob(env, &lp);
   CPXcloseCPLEX(&env);

}
