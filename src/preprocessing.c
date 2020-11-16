#include "scp.h"
#include "preprocessing.h"

/**
 * Preprocess the scp using the CPLEX library.
 * 
 * 
 * It creates the CPLEX environment and instantiate the CPLEX problem.
 * It then builds the model, and solves it at the root node, extracting the preprocessed problem after.
 * 
 * @param inst instance of the scp
 * @returns 0 if executed successfully
*/
int scpPreprocessing(instance *inst)
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
   
    // import problem from file
   CPXreadcopyprob(env, lp, inst->input_file, "LP");
   
   inst->num_cols = CPXgetnumcols(env, lp);
   inst->num_rows = CPXgetnumrows(env, lp);

   performPreprocessing(env, lp, inst);   
   // free and close CPLEX model
   CPXfreeprob(env, &lp);
   CPXcloseCPLEX(&env);
   return error;
}



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
