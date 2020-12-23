/**
 * This file contains functions needed for computing the solutions to the
 * scp using the CPLEX library.
 *
 * @author Petrucci Enrico
*/

#include "scp.h"
#include "genericCallbacks.h"
#include "constraintIntersection.h"


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

    if (inst->callback == 0)
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
 * Callback function called each time that Cplex decides to branch.
 * It evaluates if the CPLEX proposed branch is better than the one based on the constraints intersections.
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
    double end;
    instance *inst = (instance *)cbhandle; // casting of cbhandle
    int threadNum;
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &threadNum);

    // Default branching
    if (inst->branching == 0)
    {
        inst->defaultBranching[threadNum]++;
        end = second();
        inst->timeInCallback[threadNum] += end - start;
        return 0;
    }
    else // explore constraint branching
    {
        // only if CPLEX proposes a branching that generates 2 new nodes.
        if (nodecnt == 2)
        {
            int node;
            CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
            // sort based on the reduced costs
            if(node==0 && inst->sort==2)
            {
                CPXLPptr nodelp;
                CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
                start = second();

                
                double* reducedCosts = (double *)calloc(inst->num_cols, sizeof(double));
                
                
                CPXgetdj (env, nodelp, reducedCosts, 0, inst->num_cols-1);
                
                FILE *f;
                f = fopen("ReducedCosts.txt", "w");
                //fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
                fprint_array(f, reducedCosts, inst->num_cols);
                fclose(f);                    
             
                end = second();

                printf("Retrieved reduced costs in %f up to now %f \n", end - start, end-inst->startTime);

                start = second();
                computeConstraintScoresReducedCosts(inst, reducedCosts);
                end = second();
                
                printf("computed constraint score in %f up to now %f \n", end - start, end-inst->startTime);

                    
                f = fopen("ScoresConstr.txt", "w");
                fprint_array(f, inst->constraintScoresD, inst->numIntersections);
                fclose(f);        

                // f = fopen("ScoresVar.txt", "w");
                // fprint_array_int(f, inst->variableScores, inst->num_cols);
                // fclose(f);        

                start = second();
                
                // auxiliary arrays used in the merge sort, one for the intersections and one for their lengths
                int** aux= (int **)calloc(inst->numIntersections, sizeof(int*));
                int *aux1= (int *)calloc(inst->numIntersections, sizeof(int));
                double *aux2= (double *)calloc(inst->numIntersections, sizeof(double));

                
                merge_sort2(inst->repeatedNum, inst->numIntersections-1, aux, aux1, aux2, inst);
                
                
                end = second();
                printf("Sorted intersections wrt the score in %f up to now %f \n", end - start, end-inst->startTime);


                f = fopen("ScoresConstrSorted1.txt", "w");
                // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
                // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
                fprint_array(f, inst->constraintScoresD, inst->numIntersections);
                fclose(f);        

                free(aux);
                free(aux1);
                free(aux2);
                free(reducedCosts);

                
            }
            int foundConstraint = 0;
            int i = -1;
            double obj;

            double *x=inst->xs[threadNum];
            // get pseudocosts on the variables computed by CPLEX
            double *pseudocostDown = inst->varPseudocostsDown[threadNum];
            double *pseudocostUp = inst->varPseudocostsUp[threadNum];
            
            CPXgetcallbackpseudocosts(env, cbdata, wherefrom, pseudocostUp, pseudocostDown, 0, inst->num_cols - 1);

            CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);
            
            //printf("LP relaxation solved optimally it has the objective %f\n", obj);
            CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, inst->num_cols - 1);
            double Max = 0;
            
            // double sec1 = second(); 
            // inst->callbackTime1+=sec1-start;
            // 0 : lookahead after evaluating all the constraints. 
            // 1 : consider only variables > 0. 
            // 2 : consider only constraints containing best variable found by CPLEX
            // 3 : lookahead while managing each constraint
            switch (inst->constraintBranchVer)
            {
                case 0:
                {
                    Max = findBestBranchingConstraintFracVarLookAhead(&i, x, pseudocostDown, pseudocostUp, inst, threadNum);
                    break;
                }
                case 1:
                {
                    Max = findBestBranchingConstraintFracVar(&i, x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
                case 2:
                {
                    Max = findBestBranchingConstraintContainingVar(&i, indices[0], x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
                case 3:
                {
                    Max = findBestBranchingConstraintLookAhead(&i, x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
            }

            double varCostDown = (0.001 > pseudocostDown[indices[0]] * x[indices[0]] ? 0.001 : pseudocostDown[indices[0]] * x[indices[0]]);
            double varCostUp = (0.001 > pseudocostUp[indices[0]] * (1 - x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]] * (1 - x[indices[0]]));

            if (i != -1 && Max > inst->delta * (varCostDown * varCostUp))
            {
                addBranchingChilds(env, cbdata, wherefrom, obj, i, threadNum, inst);
                *useraction_p = CPX_CALLBACK_SET;
                inst->constraintBranching[threadNum]++;
            }
            else
            {
                inst->defaultBranching[threadNum]++;
            }
            end = second();
            inst->timeInCallback[threadNum] += end - start;
            return 0;
        }
        else
        {
            inst->defaultBranching[threadNum]++;
            end = second();
            inst->timeInCallback[threadNum] += end - start;
            return 0;
        }
    }
}




/**
 * For one constraint at a time it computes the fractional value of the variables
 * and the relative estimate for the score. It stops when the value does not get 
 * better after a certain number of constraints evaluated (look-ahead).
 * 
 * 
 * @param n index of the most promesing constraint found
 * @param x solution of the lp relaxation
 * @param pseudocostDown array containing pseudocost down for the variables
 * @param pseudocostUp array containing pseudocost up for the variables
 * @param inst instance of the scp
 * 
 * 
 * @returns Max the estimate of the product score
*/
double findBestBranchingConstraintLookAhead(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double estimateScoreDown;
    double estimateScoreUp;
    double constraintScorePrevision;

    double epsilon = 0.0001;

    double Max = 0;
    int best_i = -1;
    double sum = 0;
    int numIntersections= inst->numIntersections;
    int var;
    int constraintNotUsable = 0;
    int** intersections = inst->intersections;

    int nonImprovingConstraints = 0;
    int lookAhead = inst->lookAhead;
    int* intersectionsLengths = inst->intersectionsLengths; 
    int constrLength;
    int* constraint;
    double candidate;
    
    //cycle all the intersections
    for (int i = 0; i < numIntersections && nonImprovingConstraints < lookAhead; i++)
    {
        estimateScoreDown = 0;
        estimateScoreUp = INT_MAX;
        sum = 0;
        constrLength = intersectionsLengths[i];
        constraint = intersections[i];
            
        // for each variable in the current constraint update sum ad compute pseudocosts
        for (int k = 0; k < constrLength; k++)
        {
            var = constraint[k];
            
            sum += x[var];
            estimateScoreDown += pseudocostDown[var] * x[var];
            candidate = pseudocostUp[var];
            
            if(estimateScoreUp > candidate)
                    estimateScoreUp = candidate;
        }

        estimateScoreUp *= (1 - sum);

        if (sum > 1e-05 * inst->intersectionsLengths[i] && sum < 1 - 1e-05 * inst->intersectionsLengths[i]) // set covering
        {
            constraintScorePrevision = estimateScoreDown * estimateScoreUp;
            
            if (Max < constraintScorePrevision)
            {
                Max = constraintScorePrevision;
                best_i = i;
                nonImprovingConstraints=0;
            }
            else
            {
                nonImprovingConstraints++;   
            }
        }        
    }
    *n = best_i;

    return (Max);
}



/**
 * For each variable that is > 0 it updates the estimates for each of the constraints in which that variable appears.
 * It also computes in the same way the sum of the variables for each constrain in order to evaluate 
 * if they are usable for branching or not.
 * 
 *
 * @param n index of the most promesing constraint found
 * @param x solution of the lp relaxation
 * @param pseudocostDown array containing pseudocost down for the variables
 * @param pseudocostUp array containing pseudocost up for the variables
 * @param inst instance of the scp
 * 
 * 
 * @returns Max the estimate of the product score
 */
double findBestBranchingConstraintFracVar(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double start = second();
    int numIntersections = inst->numIntersections;
    
    //prepare arrays for containing the estimates and the sum for each constraint 
    double *estimateScoreDown = (double *)calloc(numIntersections, sizeof(double));
    double *estimateScoreUp = (double *)malloc(numIntersections * sizeof(double));
    // initialize estimateScoreUp to a high value in a fast way
    memset(estimateScoreUp, 127, numIntersections*sizeof(double));

    double *sum = (double *)calloc(numIntersections, sizeof(double));
    double constraintScorePrevision;
    
    double epsilon = 1e-05;

    int numVariables = inst->num_cols;
    int ** varConstrTable = inst->varConstrTable; // tells for each variable in which constraints it is
    
    int* constraintCounters = inst->constraintCounter; // array containing for each variable the number of constraints

    double end = second();
    inst->initTime+=end-start;
    start = end;
    int* constraintsWithVar; // array containing for a single variable the constraints in which it is present
    int constraint;
    double candidate; // for the estimate of the score up
    int constraintCounter;

    // cycle through all the variables
    for (int i = 0; i < numVariables; i++)
    {
        if (x[i] > epsilon)
        {
            constraintsWithVar = varConstrTable[i];
            constraintCounter = constraintCounters[i];
            // cycle all the constraints in which the variable is present
            for (int j = 0; j < constraintCounter; j++)
            {   
                constraint = constraintsWithVar[j];
                sum[constraint] += x[i];
                estimateScoreDown[constraint] += pseudocostDown[i] * x[i];
                candidate = pseudocostUp[i];
                if(estimateScoreUp[constraint] > candidate)
                    estimateScoreUp[constraint] = candidate;
            }
        }
    }

    end=second(); 
    inst->firstLoop+=end-start;
    start = end;
    double Max = 0;
    int best_i = -1;
    int* intersectionsLengths = inst->intersectionsLengths; 
    int constrLength;
    double currSum;
    double estimateScoreUpComplete;

    for (int i = 0; i < numIntersections; i++)
    {
        constrLength = intersectionsLengths[i];
        currSum = sum[i];
        // constraint is usable
        if (currSum > epsilon * constrLength && currSum < 1 - epsilon * constrLength)
        {
            // multiply the minimum pseudocode up for the fraction that is needed for the sum to reach 1
            estimateScoreUpComplete = estimateScoreUp[i] * (1 - currSum);
            constraintScorePrevision = estimateScoreUpComplete * estimateScoreDown[i];
            if (Max < constraintScorePrevision)
            {
                Max = constraintScorePrevision;
                best_i = i;
                //printf("new max for Constraints prevision: %f\n", Max);
            }
        }   
    }

    *n = best_i;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(sum);

    end=second(); 
    inst->secondLoop+=end-start;
    
    return (Max);
}


/**
 * For each variable that is > 0 it computes the sum of the variables for each constrain in order to evaluate 
 * if they are usable for branching or not.
 * Only for the constraints that are usable it computes the estimates and keep going until there is no
 * improvement after a certain number of constraints from the last best found (look-ahead)
 * 
 * 
 *
 * @param n index of the most promesing constraint found
 * @param x solution of the lp relaxation
 * @param pseudocostDown array containing pseudocost down for the variables
 * @param pseudocostUp array containing pseudocost up for the variables
 * @param inst instance of the scp
 * @param threadNum used for accessing the preallocated area for the current thread
 * 
 * 
 * @returns Max the estimate of the product score
 */
double findBestBranchingConstraintFracVarLookAhead(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst, int threadNum)
{
    //printf("in callback with %d\n", threadNum);
    double start = second();
    int numIntersections = inst->numIntersections;
    
    //prepare arrays for containing the estimates and the sum for each constraint 
    // initialize estimateScoreUp to a high value in a fast way
    double *sum = inst->sum[threadNum];
    memset(sum, 0, numIntersections*sizeof(double));
    double constraintScorePrevision;
    double epsilon = 1e-05;
    int numVariables = inst->num_cols;
    int ** varConstrTable = inst->varConstrTable; // tells for each variable in which constraints it is
    int* constraintCounters = inst->constraintCounter; // array containing for each variable the number of constraints
    double end = second();
    inst->initTime+=end-start;
    start = end;
    int* constraintsWithVar; // array containing for a single variable the constraints in which it is present
    double candidate; // for the estimate of the score up
    int constraintCounter;

    // cycle through all the variables
    for (int i = 0; i < numVariables; i++)
    {
        if (x[i] > epsilon)
        {
            constraintsWithVar = varConstrTable[i];
            constraintCounter = constraintCounters[i];
            // cycle all the constraints in which the variable is present
            for (int j = 0; j < constraintCounter; j++)
            {   
                sum[constraintsWithVar[j]] += x[i];
            }
        }
    }

    end=second(); 
    inst->firstLoop+=end-start;
    start = end;
    double Max = 0;
    int best_i = -1;
    int* intersectionsLengths = inst->intersectionsLengths; 
    int constrLength;
    double currSum;
    double estimateScoreUpComplete;
    double estimateScoreDown; 
    double estimateScoreUp;        

    int* constraint;
    
    int var;
    int constraintNotUsable = 0;
    int** intersections = inst->intersections;
    int nonImprovingConstraints = 0;
    int lookAhead = inst->lookAhead;
    //cycle all the intersections
    
    for (int i = 0; i < numIntersections && nonImprovingConstraints < lookAhead; i++)
    {
        estimateScoreDown = 0;
        estimateScoreUp = INT_MAX;
        
        constrLength = intersectionsLengths[i];
        currSum = sum[i];
         
        // constraint is usable
        if (currSum > epsilon * constrLength && currSum < 1 - epsilon * constrLength)
        {
            constraint = intersections[i];
            // for each variable in the current constraint update sum ad compute pseudocosts
            for (int k = 0; k < constrLength; k++)
            {
                var = constraint[k];
                estimateScoreDown += pseudocostDown[var] * x[var];

                candidate = pseudocostUp[var];

                if(estimateScoreUp > candidate)
                    estimateScoreUp = candidate;
            }

            // multiply the minimum pseudocode up for the fraction that is needed for the sum to reach 1
            estimateScoreUpComplete = estimateScoreUp * (1 - currSum);
            
            constraintScorePrevision = estimateScoreUpComplete * estimateScoreDown;
            if (Max < constraintScorePrevision)
            {
                Max = constraintScorePrevision;
                best_i = i;
                nonImprovingConstraints=0;

                //printf("new max for Constraints prevision: %f\n", Max);
            }
            else
            {
                nonImprovingConstraints++;   
            }
        }
        else
        {
            constraintNotUsable++;   
        }
    }
    // printf("constraintNotUsable %d out of %d\n",constraintNotUsable, inst->numIntersections);
    *n = best_i;
    end=second(); 
    inst->secondLoop+=end-start;
    return (Max);
}


/**
 * For each branching constraints that contains the best branching variable identified by CPLEX
 * it computes the estimates of the product scores, starting from the pseudocosts for each variable
 * which are provided by CPLEX 
 *
 *
 * @param n index of the most promesing constraint found
 * @param bestvar index of the variable chosen by CPLEX for branching
 * @param x solution of the lp relaxation
 * @param pseudocostDown array containing pseudocost down for the variables
 * @param pseudocostUp array containing pseudocost up for the variables
 * @param inst instance of the scp
 * 
 * 
 * @returns Max the estimate of the product score
 */

double findBestBranchingConstraintContainingVar(int *n, int bestVar, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double constraintScorePrevision;
    
    double epsilon = 1e-05;

    double Max = 0;
    int best_i = -1;

    int *constraintCounter = inst->constraintCounter; 
    int **varConstrTable = inst->varConstrTable;
    int **intersections=inst->intersections;
    int *intersectionsLengths = inst->intersectionsLengths;
    
    double estimateDown;
    double estimateUp;
    int constCounter = constraintCounter[bestVar]; 

    // cycle the constraints that contain the best variable
    for (int i = 0; i < constCounter; i++)
    {
        constraintScorePrevision=0;
        
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        estimateDown = 0;
        estimateUp = INT_MAX;
        
        int constraintIndex = varConstrTable[bestVar][i];
        int *constraint = intersections[constraintIndex];    
        int constraintLength = intersectionsLengths[constraintIndex];

        for (int k = 0; k < constraintLength; k++)
        {
            int currVar = constraint[k];
            sum += x[currVar];
            estimateDown += pseudocostDown[currVar] * x[currVar];
            // printf("current downpseudocost = %f\n", pseudoDown);
            double candidate = pseudocostUp[currVar];
            if(estimateUp > candidate)
                estimateUp = candidate;
        }

        estimateUp *= (1 - sum);

        if (sum > epsilon * constraintLength && sum < 1 - epsilon * constraintLength) // set covering
        {
            constraintScorePrevision = estimateDown * estimateUp;
        }
        if (Max < constraintScorePrevision)
        {
            Max = constraintScorePrevision;
            best_i = constraintIndex;
        }
    }
    *n = best_i;
    return (Max);
}


/** while using LegacyCallback manages the case in which the branching is done by using CPLEX default or the branching constraints.
 * In the second case some arrays are populated in order to contain the data required during the callbacks.
 * 
 * @param env CPLEX enviroment
 * @param lp CPLEX problem 
 * @param inst instance of the scp
 * 
 */
void solveUsingLegacyCallback(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    char str[80];

    // prepare logfile
    if (inst->branching == 0)
        sprintf(str, "logfile_defaultBranchinglegacy.txt");
    else
        sprintf(str, "logfile_constraintBranchinglegacy.txt");
    CPXsetlogfilename(env, str, "w");


    // set callback
    int plotTree=0;
    if(plotTree)
    {
        resetPlotTree(inst);
        CPXsetbranchcallbackfunc(env, plotTreeCallback, inst);
    }
    else
        CPXsetbranchcallbackfunc(env, legacyBranchingCallback, inst);
    if (inst->threads == 0)
        CPXgetnumcores(env, &inst->threads);
    CPXsetintparam(env, CPX_PARAM_THREADS, inst->threads);
    printf("Using %d threads\n", inst->threads);
    
    // during callbacks do not use reduced representation
    CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
    
    inst->constraintBranching = (int *)calloc(inst->threads, sizeof(int));
    inst->defaultBranching = (int *)calloc(inst->threads, sizeof(int));

    inst->timeInCallback = (double *)calloc(inst->threads, sizeof(double));
    // inst->timeFindingConstraint = (double *)calloc(inst->threads, sizeof(double));

    if (inst->branching == 1)
    {

        //CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_PSEUDO);
        // Get the rows of the problem
        int nnz;
        int *izero = (int *)malloc(inst->num_rows * sizeof(int));
        int *indexes = (int *)malloc(inst->num_rows * inst->num_cols * sizeof(int));         
        double *values = (double *)malloc(inst->num_rows * inst->num_cols * sizeof(double)); 
        int surplus_p;

        CPXgetrows(env, lp, &nnz, izero, indexes, values, inst->num_rows * inst->num_cols, &surplus_p, 0, inst->num_rows - 1);

        printf("pre preprocessing: %f\n", second()-inst->startTime);
        double start;
        double end;

        start = second();
        populateIntersectionsOf2(izero, indexes, nnz, inst);
        end = second();
        printf("Found %d constraint intersections non reordered in %f. Up to now %f\n", inst->numIntersections, end - start, end-inst->startTime);
        
        // start = second();
        // populateIntersectionsOf2Original(izero, indexes, nnz, inst);
        // end = second();
        // printf("Found %d constraint intersections non reordered original in %f. Up to now %f\n", inst->numIntersections, end - start, end-inst->startTime);
        
  
        start = second();
        
        // auxiliary arrays used in the merge sort, one for the intersections and one for their lengths
        int **aux= (int **)calloc(inst->numIntersections, sizeof(int*));
        int *aux1= (int *)calloc(inst->numIntersections, sizeof(int));
        merge_sort(0, inst->numIntersections-1, aux, aux1, inst);
        
        free(aux);
        free(aux1);

        end = second();
        printf("Sorted intersections in %f up to now %f \n", end - start, end-inst->startTime);
        
        printf("repeatedFirst %d\n", inst->repeatedFirst);
        if(inst->repeatedFirst==1)
        {
            start = second();
            purgeDuplicatesRepeatedFirst(inst);
            end = second();
            printf("Eliminated duplicates in %f. Repeated brougth at the top. Remaining %d constraints up to now %f \n", end - start, inst->numIntersections, end-inst->startTime);
        }
        else
        {
            start = second();
            purgeDuplicates(inst);
            end = second();
            printf("Eliminated duplicates in %f. Remaining %d constraints up to now %f \n", end - start, inst->numIntersections, end-inst->startTime);
        }
        
        start = second();
        populateVariableConstraintTable(inst);
        end = second();
        printf("Populated variable-constraint table in %f up to now %f \n", end - start, end-inst->startTime);

        if(inst->sort==1)
        {
            start = second();
            computeVariableFrequencies(indexes, nnz, inst);
            end = second();
            printf("computed variable frequencies in %f up to now %f \n", end - start, end-inst->startTime);

            start = second();
            computeConstraintScoresFreq(inst);
            end = second();
            printf("computed constraint score in %f up to now %f \n", end - start, end-inst->startTime);

                
            // FILE *f;
            // f = fopen("ScoresConstr.txt", "w");
            // fprint_array_int(f, inst->constraintScores, inst->numIntersections);
            // fclose(f);        

            // f = fopen("ScoresVar.txt", "w");
            // fprint_array_int(f, inst->variableScores, inst->num_cols);
            // fclose(f);        

            start = second();
            
            // auxiliary arrays used in the merge sort, one for the intersections and one for their lengths
            aux= (int **)calloc(inst->numIntersections, sizeof(int*));
            aux1= (int *)calloc(inst->numIntersections, sizeof(int));
            int *aux2= (int *)calloc(inst->numIntersections, sizeof(int));

            merge_sort1(inst->repeatedNum, inst->numIntersections-1, aux, aux1, aux2, inst);
            
            
            end = second();
            printf("Sorted intersections wrt the score in %f up to now %f \n", end - start, end-inst->startTime);


            // f = fopen("ScoresConstrSorted1.txt", "w");
            // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
            // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
            // fprint_array_int(f, inst->constraintScores, inst->numIntersections);
            // fclose(f);        

            free(aux);
            free(aux1);
            free(aux2);
        }
        // if(inst->sort==2)
        // {
        //     start = second();
        //     CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, 0);
        //     if (CPXmipopt(env, lp)) 
        //         print_error(" Problems on CPXmipopt");
            
        //     end = second();
        //     printf("Root node solved in %f up to now %f\n", end - start, end-inst->startTime);
            
        //     double* reducedCosts = (double *)calloc(inst->num_cols, sizeof(double));
            
        //     CPXgetdj (env, lp, reducedCosts, 0, inst->num_cols-1);
            
        //     FILE *f;
            
        //     f = fopen("ReducedCosts.txt", "w");
        //     //fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
        //     fprint_array(f, reducedCosts, inst->num_cols);

        //     fclose(f);                    
        //     CPXsetlongparam(env, CPXPARAM_MIP_Limits_Nodes, 9223372036800000000);
            

        // }

        
            // FILE *f;
            // f = fopen("ScoresConstrSorted1.txt", "w");
            // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
            // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
            // fclose(f);        

        // allocate arrays only once for each thread to use inside the callback

        inst->xs = (double **)malloc(inst->threads * sizeof(double*));
        inst->varPseudocostsDown = (double **)malloc(inst->threads * sizeof(double*));
        inst->varPseudocostsUp = (double **)malloc(inst->threads * sizeof(double*));
        inst->sum = (double **)malloc(inst->threads * sizeof(double*));
    
        for (int i=0; i<inst->threads; i++)
        {
            inst->xs[i] = (double *)malloc(inst->num_cols * sizeof(double));
            inst->varPseudocostsDown[i] = (double *)malloc(inst->num_cols * sizeof(double));
            inst->varPseudocostsUp[i] = (double *)malloc(inst->num_cols * sizeof(double));
            inst->sum[i] = (double *)malloc(inst->numIntersections * sizeof(double));
        }

        inst->values =(double *) malloc(inst->maxConstrLen * sizeof(double));
        inst->lu = (char *) malloc(inst->maxConstrLen * sizeof(char));
        inst->bd =(double *) calloc(inst->maxConstrLen, sizeof(double));
        
        for (int n = 0; n < inst->maxConstrLen; n++)
        {
            inst->values[n] = 1;
            inst->lu[n] = 'U';

        }

        free(izero);
        free(indexes);
        free(values);
    }

    inst->prepTime = second() - inst->startTime;

    // solve the problem
    if (CPXmipopt(env, lp)) 
        print_error(" Problems on CPXmipopt");
    
    inst->solution = (double *)calloc(inst->num_cols, sizeof(double));
    double obj;

    if (CPXgetx(env, lp, inst->solution, 0, inst->num_cols - 1) == 0 && CPXgetbestobjval(env, lp, &obj) == 0)
        printf("Solution found with objval = %f\n", obj);
    else
        print_error("Problems in getting the solution");

    inst->endTime=second();
    inst->executionTime = (inst->endTime - inst->startTime);
        
    inst->totalConstraintBranching = 0;
    inst->totalVariableBranching = 0;

    double callbackTime = 0;
    double findingConstraintTime = 0;

    for (int i = 0; i < inst->threads; i++)
    {
        inst->totalConstraintBranching += (inst->constraintBranching[i]);
        inst->totalVariableBranching += (inst->defaultBranching[i]);
        callbackTime += (inst->timeInCallback[i]);
    }
    if(inst->storeResults)
    {
        CPXgetcutoff( env, lp, &(inst->bestInt));
        inst->bestVal = obj;
        CPXgetmiprelgap (env, lp, &(inst->MIPgap));
        inst->exploredNodes = CPXgetnodecnt (env, lp);
        inst->remainingNodes = CPXgetnodeleftcnt (env, lp);

        saveComputationResults(inst);
    }
    
    printf("%d constraint branching\n%d default branching\n", inst->totalConstraintBranching, inst->totalVariableBranching);
    
    printf("Total time spent during preprocessing = %f which is %f%% of the total\n",inst->prepTime, 100 * inst->prepTime / inst->executionTime);
    
    printf("Total time spent in callback = %f which is %f%% of the total\n", callbackTime, 100 * callbackTime / inst->executionTime);
    // printf("Time spent in callback choosing the best constraint: %f\n", findingConstraintTime);
    // if(inst->constraintBranchVer==0)
    // {
    //     printf("Total time spent during initialization = %f which is %f%% of the total\n",inst->initTime, 100 * inst->initTime / (second() - inst->startTime));
    //     printf("Total time spent in the first loop = %f which is %f%% of the total\n",inst->firstLoop, 100 * inst->firstLoop / (second() - inst->startTime));
    //     printf("Total time spent in the second loop = %f which is %f%% of the total\n",inst->secondLoop, 100 * inst->secondLoop / (second() - inst->startTime));
    // }
    // printf("Total time spent callback start %f which is %f%% of the total\n", inst->callbackTime1, 100 * inst->callbackTime1 / (second() - inst->startTime));
    // printf("Total time spent finding the constraint for branching = %f which is %f%% of the total\n", inst->callbackTime2, 100 * inst->callbackTime2 / (second() - inst->startTime));
    // printf("Total time spent in building the childs = %f which is %f%% of the total\n", inst->callbackTime3, 100 * inst->callbackTime3 / (second() - inst->startTime));
    

    if (inst->branching == 1)
    {
        for (int i = 0; i < inst->numIntersections; i++)
        {
            free(inst->intersections[i]);
            //printf("Free intersection %d of %d\n", i, inst->numIntersections);
            
        }
        
        printf("Freed all the intersections\n");
        
        free(inst->intersections);
        free(inst->intersectionsLengths);
        printf("Free intersections\n");
        for(int i = 0; i < inst->num_cols; i++)
        {
            free(inst->varConstrTable[i]);
        }
        free(inst->varConstrTable);
        free(inst->constraintCounter);
        printf("Free varConstrTable\n");
        for (int i=0; i<inst->threads; i++)
        {
            free(inst->xs[i]);
            free(inst->varPseudocostsDown[i]);
            free(inst->varPseudocostsUp[i]);
            free(inst->sum[i]);
        }
        free(inst->xs);
        free(inst->varPseudocostsDown);
        free(inst->varPseudocostsUp);
        free(inst->sum);
        free(inst->values);
        free(inst->bd);
        free(inst->lu);
        if (inst->sort==1)
        {
            free(inst->constraintScores);
            free(inst->variableScores);
        }
        if(inst->sort==2)
        {
            free(inst->constraintScoresD);
        }
    }
    
    free(inst->solution);
    free(inst->constraintBranching);
    free(inst->defaultBranching);
    free(inst->timeInCallback);
    // free(inst->timeFindingConstraint);
}


/** while using LegacyCallback manages the case in which the branching is done by using CPLEX default or the branching constraints.
 * In the second case some arrays are populated in order to contain the data required during the callbacks.
 * 
 * @param env CPLEX enviroment
 * @param cbdata info for the callback
 * @param wherefrom info for the callback
 * @param obj objective value
 * @param i constraint to use
 * @param inst instance of the scp
 * 
 */
void addBranchingChilds(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, int threadNum, instance *inst)
{
    // generating information for adding the constraint to CPLEX
    // in the first case the sum of the variables is set equal to zero,
    // in the second case the sum of the variables is set to be equal or greater than 1
    double rhs1 = 0;
    double rhs2 = 1;
    char sense1 = 'E';
    char sense2 = 'G';
    int izero = 0;
    // print_array_int(inst->intersections[i], inst->intersectionsLengths[i]);
    // print_array(inst->bd, inst->intersectionsLengths[i]);
    // print_array_char(inst->lu, inst->intersectionsLengths[i]);
    /* We want to branch. */
    int child1, child2;
    CPXbranchcallbackbranchbds(env, cbdata, wherefrom, inst->intersectionsLengths[i], inst->intersections[i], inst->lu, inst->bd, obj, NULL, &child1);
    //CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->intersectionsLengths[i], &rhs1, "E", &izero, inst->intersections[i], inst->values, obj, NULL, &child1);
    CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->intersectionsLengths[i], &rhs2, "G", &izero, inst->intersections[i], inst->values, obj, NULL, &child2);
}

/**
 * Save the details of the computations in a .csv file
 * @param inst instance of the scp
 * 
 */

void saveComputationResults(instance *inst)
{
    char fileName[1000];
    char folder[1000];
    int length = strlen(inst->input_file);
    int start = 21;
    int end = length-3;
    char input_file[100];
    int i=0;
    for(; i<end-start; i++)
    {
        input_file[i]=inst->input_file[i+start];
    }
    input_file[i]=0;
    sprintf(folder, "../results/%d_%d_%d_%.1f_%d_%d_%d_%d_%d", inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->reverse, inst->repeatedFirst, inst->sort, inst->average);
    sprintf(fileName, "../results/%d_%d_%d_%.1f_%d_%d_%d_%d_%d/%s_%d_%d_%d_%.1f_%d_%d_%.0f_%d_%d_%d_%d.csv", inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->reverse, inst->repeatedFirst, inst->sort, inst->average, input_file, inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->threads, inst->timelimit, inst->reverse, inst->repeatedFirst, inst->sort, inst->average);
    
    
    printf("Saving results in file %s\n", fileName);
    // create folder if not present
    struct stat st = {0};
    if (stat(folder, &st) == -1) {
        mkdir(folder, 0700);
    }

    FILE *f;
    if( access( fileName, F_OK ) != -1 ) // file exists
    {
        f = fopen(fileName, "a"); 
    }
    else // file doesn't exist
    {
        f = fopen(fileName, "w");
        fprintf(f, "Instance,Time,Best Int.,Best Val.,MIP Gap,Nodes,Nodes Left,Constraint Branching,Varaible Branching,Seed,Threads,Callback,Branching,constraintBranchVer,delta,lookAhead,Reverse,RepeatedFirst,Sort,Average,PreprocessingPercentage\n"); 
    }
    fprintf(f, "%s,%f,%f,%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%d,%.1f,%d,%d,%d,%d,%d,%f\n",inst->input_file, inst->executionTime, inst->bestInt, inst->bestVal, inst->MIPgap, inst->exploredNodes, inst->remainingNodes, inst->totalConstraintBranching, inst->totalVariableBranching, inst->seed, inst->threads, inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->reverse, inst->repeatedFirst, inst->sort, inst->average,100 * inst->prepTime / inst->executionTime); 
    fclose(f);
}
    

/**
 * For each possible branching constraints it computes the estimates of the product scores,
 * starting from the pseudocosts for each variable which are provided by CPLEX 
 *
 *
 * @param n index of the most promesing constraint found
 * @param x solution of the lp relaxation
 * @param pseudocostDown array containing pseudocost down for the variables
 * @param pseudocostUp array containing pseudocost up for the variables
 * @param inst instance of the scp
 * 
 * 
 * @returns Max the estimate of the product score
*/
double findBestBranchingConstraint(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double *estimateScoreDown = (double *)calloc(inst->numIntersections, sizeof(double));
    double *estimateScoreUp = (double *)calloc(inst->numIntersections, sizeof(double));
    double *constraintScorePrevision = (double *)calloc(inst->numIntersections, sizeof(double));
    double epsilon = 0.0001;

    double Max = 0;
    int best_i = -1;
    
    // cycle the intersections
    for (int i = 0; i < inst->numIntersections; i++)
    {
        double pseudoDown = 0;
        double pseudoUp = INT_MAX;
        double sum = 0;
        
        // for each variable in the current constraint update sum ad compute pseudocosts
        for (int k = 0; k < inst->intersectionsLengths[i]; k++)
        {
            sum += x[inst->intersections[i][k]];
            pseudoDown += pseudocostDown[inst->intersections[i][k]] * x[inst->intersections[i][k]];

            double candidate = pseudocostUp[inst->intersections[i][k]];
            pseudoUp = (pseudoUp > candidate ? candidate : pseudoUp);
        }

        pseudoUp *= (1 - sum);

        if (sum > 1e-05 * inst->intersectionsLengths[i] && sum < 1 - 1e-05 * inst->intersectionsLengths[i]) // set covering
        {
            // printf("breakpoint");
            estimateScoreDown[i] = pseudoDown;
            estimateScoreUp[i] = pseudoUp;
            constraintScorePrevision[i] = pseudoDown * pseudoUp;
        }
        else
        {
            estimateScoreDown[i] = epsilon;
            estimateScoreUp[i] = epsilon;
            constraintScorePrevision[i] = epsilon * epsilon;
        }
        if (Max < constraintScorePrevision[i])
        {
            Max = constraintScorePrevision[i];
            best_i = i;
        }
    }
    *n = best_i;

    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return (Max);
}



/** 
 * Attempt for plotting the tree structure of the decision tree
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
int plotTreeCallback(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                            const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
{
    double start = second();
    instance *inst = (instance *)cbhandle; // casting of cbhandle
    int threadNum;
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &threadNum);
    int node;
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
    printf("node = %d\n", node);

    // int indices0[bdcnt];
    // int indices1[bdcnt];

    // for(int i=0; i<bdcnt; i++)
    // {
    //     if(i<nodebeg[1])
    //     {
            
    //     }
    //     else
    //     {

    //     }

    // }
    // printf("bdcnt = %d\n", bdcnt);
    // print_array_int(nodebeg, nodecnt);
    // print_array_char(lu, bdcnt);
    int seqnum0;
    int seqnum1;
    if(nodecnt==2)
    {
        CPXbranchcallbackbranchbds( env, cbdata, wherefrom, nodebeg[1], &(indices[nodebeg[0]]), &(lu[nodebeg[0]]), &(bd[nodebeg[0]]), nodeest[0], NULL, &seqnum0);
        CPXbranchcallbackbranchbds( env, cbdata, wherefrom, bdcnt-nodebeg[1], &(indices[nodebeg[1]]), &(lu[nodebeg[1]]), &(bd[nodebeg[1]]) , nodeest[1], NULL, &seqnum1);
        printf("node0 = %d, node1 = %d\n", seqnum0,seqnum1);
        *useraction_p = CPX_CALLBACK_SET;
    
        char str[1000];
        int length = strlen(inst->input_file);
        int start = 21;
        int end = length-3;
        char input_file[100];
        int i=0;

        for(; i<end-start; i++)
        {
            input_file[i]=inst->input_file[i+start];
        }

        input_file[i]=0;
        sprintf(str, "../TreePlots/tree%s.txt", input_file);
        FILE *f;
        f = fopen(str, "a");
        fprintf(f,"%d %d \n",node, seqnum0);
        
        fprintf(f,"%d %d \n",node, seqnum1);
        fclose(f);
    }
    return 0;
}



void resetPlotTree(instance* inst)
{
    char str[1000];
	int length = strlen(inst->input_file);
	int start = 21;
	int end = length-3;
	char input_file[100];
	int i=0;

	for(; i<end-start; i++)
	{
		input_file[i]=inst->input_file[i+start];
	}

	input_file[i]=0;
	sprintf(str, "../TreePlots/tree%s.txt", input_file);
	FILE *f;
	f = fopen(str, "w");
	fclose(f);
}