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
    
    // int depth;
    // CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
    // printf("depth = %d. ", depth);
    
    // Default branching
    if (inst->branching == 0)
    {
        inst->defaultBranching[threadNum]++;
        // printf("default branching\n");
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
                sortIntersectionWRTreducedCosts(env, cbdata, wherefrom, inst);
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
            // 4 : strong branching on constraints containing CPLEX proposed variable.
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
                case 4:
                {
                    Max = findBestBranchingConstraintStrongContVar(env, cbdata, wherefrom, &i, indices[0], threadNum, obj, inst);
                    break;
                }
                case 5:
                {
                    Max = findBestBranchingConstraintStrongContVar(env, cbdata, wherefrom, &i, indices[0], threadNum, obj, inst);
                    break;
                }
            }
            
            double varProductScore;

            if(inst->constraintBranchVer >= 4)
            {
                varProductScore = computeVariableProductScore(inst, x, obj, threadNum, bdcnt, nodebeg, indices, lu, bd);
            }
            else
            {
                double varCostDown = (0.001 > pseudocostDown[indices[0]] * x[indices[0]] ? 0.001 : pseudocostDown[indices[0]] * x[indices[0]]);
                double varCostUp = (0.001 > pseudocostUp[indices[0]] * (1 - x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]] * (1 - x[indices[0]]));
                varProductScore = (varCostDown * varCostUp);
            }
            
            if (i != -1 && Max > inst->delta * (varProductScore))
            {
                addBranchingChilds(env, cbdata, wherefrom, obj, i, threadNum, inst);
                *useraction_p = CPX_CALLBACK_SET;
                inst->constraintBranching[threadNum]++;
                // printf("Constraint branching\n");
                // inst->constraintBranchingDepth[depth]++;
                
            }
            else
            {
                inst->defaultBranching[threadNum]++;
                // printf("default branching\n");
                // inst->variableBranchingDepth[depth]++;
                
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
            // printf("default branching\n");
            // inst->variableBranchingDepth[depth]++;
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


/**
 * For each branching constraints that contains the best branching variable identified by CPLEX
 * it computes the strong branching, returns the product score of the best constraint and the index of the constraint. 
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

double findBestBranchingConstraintStrongContVar(CPXCENVptr env, void *cbdata ,int wherefrom, int *n, int bestVar, int threadNum, double obj, instance *inst)
{
    //printf("best variable: %d\n", bestVar);

    double epsilon = 1e-05;
    double Max = 0;
    int best_i = -1;
    
    double currScore;

    double* x = inst->xs[threadNum]; 

    //clone lp and prepare for solving lp relaxations

    // strong branching
    CPXLPptr nodelp;
    CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);
    
    int cur_numcols = CPXgetnumcols (env, nodelp);
    int cur_numrows = CPXgetnumrows (env, nodelp);
    
    getBase(env, nodelp, inst, threadNum, cur_numcols, cur_numrows);          
            
    int* cstat = inst->cstats[threadNum];
    int* rstat = inst->rstats[threadNum];
    int status;
    
    inst->lps[threadNum] = CPXcloneprob (inst->envs[threadNum], nodelp, &status);                
    // printf("clone status %d\n", status);

    // state that the problem is a lp problem
    CPXchgprobtype(inst->envs[threadNum], inst->lps[threadNum], 0);
    
    CPXgetub(inst->envs[threadNum], inst->lps[threadNum], inst->ogBd[threadNum], 0, cur_numcols-1);

    
    int *constraintCounter = inst->constraintCounter; 
    int **varConstrTable = inst->varConstrTable;
    int **intersections = inst->intersections;
    int *intersectionsLengths = inst->intersectionsLengths;
    
    int constCounter = constraintCounter[bestVar]; 
    //printf("Number of constraint containing the variable = %d\n",constCounter);
    // cycle the constraints that contain the best variable
    int useConstraintsContainingVar = 1;
    if(useConstraintsContainingVar == 1)
    {
        for (int i = 0; i < constCounter; i++)
        {
            double sum = 0;
            // cycle all the variables in the intersection and get the sum
            
            int constraintIndex = varConstrTable[bestVar][i];
            int *constraint = intersections[constraintIndex];
            int constraintLength = intersectionsLengths[constraintIndex];
            
            // printf("considering constraint :");
            // print_array_int(constraint, constraintLength);
            
            for (int k = 0; k < constraintLength; k++)
            {
                sum += x[constraint[k]];
            }
            
            if (sum > epsilon * constraintLength && sum < 1 - epsilon * constraintLength) // set covering
            {
                currScore = computeStrongBranchingOnConstraint(inst, obj, cstat, rstat, cur_numcols, threadNum, constraintIndex);
                //printf("currScore = %f\n", currScore);
                if (Max < currScore)
                {
                    Max = currScore;
                    best_i = constraintIndex;
                }
            }       
        }
    }
    else
    {
        for (int i = 0; i < inst->numIntersections; i++)
        {
            double sum = 0;
            // cycle all the variables in the intersection and get the sum
            
            int *constraint = intersections[i];
            int constraintLength = intersectionsLengths[i];
            
            // printf("considering constraint :");
            // print_array_int(constraint, constraintLength);
            
            for (int k = 0; k < constraintLength; k++)
            {
                sum += x[constraint[k]];
            }
            
            if (sum > epsilon * constraintLength && sum < 1 - epsilon * constraintLength) // set covering
            {
                currScore = computeStrongBranchingOnConstraint(inst, obj, cstat, rstat, cur_numcols, threadNum, i);
                //printf("currScore = %f\n", currScore);
                if (Max < currScore)
                {
                    Max = currScore;
                    best_i = i;
                }
            }       
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
    switch(inst->branching)
    {
        case 0:
            sprintf(str, "logfile_defaultBranchinglegacy.txt");
            break;
        case 1:
            sprintf(str, "logfile_ConstraintBranchinglegacy%d.txt", inst->constraintBranchVer);
            break;
        case 2:
            sprintf(str, "logfile_defaultStrongBranchinglegacy.txt");
            break;
        case 3:
            sprintf(str, "logfile_constraintStrongBranchinglegacy.txt");
            break;
    }
    
    CPXsetlogfilename(env, str, "w");

    if(inst->branching==0 || inst->branching==1 )
    {
        CPXsetbranchcallbackfunc(env, legacyBranchingCallback, inst);
        //CPXsetintparam (env, CPXPARAM_Preprocessing_Presolve, CPX_OFF);
        
        if(inst->constraintBranchVer == 5)
        {
            CPXsetintparam (env, CPXPARAM_MIP_Limits_StrongCand, CPX_BIGINT);
            CPXsetlongparam (env, CPXPARAM_MIP_Limits_StrongIt, CPX_BIGLONG);
            CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
        }
    }
    else
    {
        CPXsetbranchcallbackfunc(env, legacyBranchingCallbackStrong, inst);
        CPXsetintparam (env, CPXPARAM_MIP_Limits_StrongCand, CPX_BIGINT);
        CPXsetlongparam (env, CPXPARAM_MIP_Limits_StrongIt, CPX_BIGLONG);
        CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);
    }
    
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

    if (inst->branching == 1 || inst->branching == 3)
        prepareBranchingConstraints(env, lp, inst);
    
  
    // if(inst->branching == 1)
    // {
    //     inst->constraintBranchingDepth = (int*)calloc(100, sizeof(int));
    //     inst->variableBranchingDepth = (int*)calloc(100, sizeof(int));
    // }
  
    if (inst->branching == 3 || inst->branching == 2|| (inst->branching == 1 && inst->constraintBranchVer >= 4))
    {

        int error;
        inst->envs = (CPXENVptr*)malloc(inst->threads * sizeof(CPXENVptr*));
        inst->lps = (CPXLPptr*)malloc(inst->threads * sizeof(CPXLPptr*));
        inst->cstats = (int**)malloc(inst->threads * sizeof(int*));
        inst->rstats = (int**)malloc(inst->threads * sizeof(int*));
        inst->cstatDims = (int*)malloc(inst->threads * sizeof(int));
        inst->rstatDims = (int*)malloc(inst->threads * sizeof(int));
        
        inst->ogBd = (double**)malloc(inst->threads * sizeof(double*));
        inst->ogIndices = (int**)malloc(inst->threads * sizeof(int*));
        inst->ogLu = (char**)malloc(inst->threads * sizeof(char*));
        
        
        for (int i=0; i<inst->threads; i++)
        {
            inst->envs[i] = CPXopenCPLEX(&error);
            printf("opened cplex enviroment for solving lp relaxation: %d\n", error);
            inst->cstatDims[i] = inst->num_cols * 2;
            inst->rstatDims[i] = inst->num_rows * 2;
            inst->cstats[i] = (int*)malloc(inst->num_cols * inst->cstatDims[i] * sizeof(int));
            inst->rstats[i] = (int*)malloc(inst->num_rows * inst->rstatDims[i] * sizeof(int));
            
            inst->ogBd[i] = (double *) malloc (inst->cstatDims[i]*sizeof(double));
            inst->ogIndices[i] = (int *) malloc (inst->cstatDims[i]*sizeof(int));
            inst->ogLu[i] = (char *) malloc (inst->cstatDims[i]*sizeof(char));

            for(int k=0; k<inst->cstatDims[i]; k++)
            {
                inst->ogIndices[i][k] = k;
                inst->ogLu[i][k] = 'U';
            } 
        }

        inst->nInfo = (nodeInfo**)malloc(inst->threads * sizeof(nodeInfo*));
	    inst->nInfoLengths = (int* )malloc(inst->threads * sizeof(int));
	    inst->nInfoIndex = (int* )calloc(inst->threads, sizeof(int));
        for (int i=0; i<inst->threads; i++)
        {
            inst->nInfoLengths[i] = 1000;
            inst->nInfo[i] = (nodeInfo*) malloc(inst->nInfoLengths[i] * sizeof(nodeInfo));
        }
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
        printf("Explored nodes = %d\n", inst->exploredNodes);
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
    // if(inst->branching == 1)
    // {
    //     FILE *f;
    //     f = fopen("branchingConstraintDepth.txt", "w");
    //     fprintf(f, "Variable branching\n");
    //     fprint_array_int(f, inst->variableBranchingDepth, 100);
        
    //     fprintf(f, "Constraint branching\n");
    //     fprint_array_int(f, inst->constraintBranchingDepth, 100);
    // }    
    
    if (inst->branching == 3 || inst->branching == 2 || (inst->branching == 1 && inst->constraintBranchVer >= 4))
    {

        //saveNodeInfoToFile(inst);
        for(int i=0; i<inst->threads; i++)
        {
            free(inst->nInfo[i]);
        }
        free(inst->nInfo);
        free(inst->nInfoLengths);
        free(inst->nInfoIndex);
  
        for (int i=0; i<inst->threads; i++)
        {
            free(inst->ogBd[i]);
            free(inst->ogIndices[i]);
            free(inst->ogLu[i]);
            free(inst->cstats[i]);
            free(inst->rstats[i]);        
            CPXcloseCPLEX(&inst->envs[i]);
        }
        free(inst->ogBd);
        free(inst->ogIndices);
        free(inst->ogLu);
        free(inst->cstats);
        free(inst->rstats);
        free(inst->lps);
        free(inst->envs);
        free(inst->cstatDims);
        free(inst->rstatDims);

        
        
    }

    if (inst->branching == 1 || inst->branching == 3)
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

/**
 * From the original formulation of the problem computes the branching constraints 
 * and in case sorts them, removes the duplicates and eventually brings the one that 
 * were duplicates on top  
 *
 * @param env Cplex environment
 * @param lp Cplex problem 
 * @param inst instance of the scp
 * 
*/

void prepareBranchingConstraints(CPXCENVptr env, CPXLPptr lp, instance *inst)
{
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

            
        start = second();
        
        // auxiliary arrays used in the merge sort, one for the intersections and one for their lengths
        aux= (int **)calloc(inst->numIntersections, sizeof(int*));
        aux1= (int *)calloc(inst->numIntersections, sizeof(int));
        int *aux2= (int *)calloc(inst->numIntersections, sizeof(int));

        merge_sort1(inst->repeatedNum, inst->numIntersections-1, aux, aux1, aux2, inst);
        
        
        end = second();
        printf("Sorted intersections wrt the score in %f up to now %f \n", end - start, end-inst->startTime);


        free(aux);
        free(aux1);
        free(aux2);
    }

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

double computeStrongBranchingOnConstraint(instance* inst, double obj, int* cstat, int* rstat, int ogNumCols, int threadNum, int i)
{
    int verbose = 0;
    
    CPXENVptr env = inst->envs[threadNum];
    CPXLPptr lp = inst->lps[threadNum]; 
    
    if(verbose)
    {
        printf("Branching on constraint:");
        print_array_int(inst->intersections[i], inst->intersectionsLengths[i]);
        CPXwriteprob(env, lp, "Model0.lp", NULL);
    }

    double epsilon = 0.001;
   
    double leftobj = 0;
    double rightobj = 0;

    double delta1;
    double delta2;
    
    double rhs2 = 1;
    char* name2 = "up";

    int izero = 0;
    
    // test 0 branch
    CPXcopybase (env, lp, cstat, rstat);
    CPXchgbds(env, lp, inst->intersectionsLengths[i], inst->intersections[i], inst->lu, inst->bd);

    if(verbose)
    {
        CPXwriteprob(env, lp, "Model1.lp", NULL);
    }

    if (CPXlpopt(env, lp)) // solve the problem
        print_error(" Problems on CPXlpopt1");

    CPXgetobjval(env, lp, &leftobj);
            
    delta1 = epsilon > leftobj-obj ? epsilon : leftobj-obj;


    //reset problem
    CPXchgbds(env, lp, ogNumCols, inst->ogIndices[threadNum], inst->ogLu[threadNum], inst->ogBd[threadNum]);
    
    if(verbose)
    {
        CPXwriteprob(env, lp, "Model2.lp", NULL);
    }

    //CPXdelrows(env, lp, branchingConstraint, branchingConstraint);

    // test 1 branch 
    CPXcopybase (env, lp, cstat, rstat);

    CPXaddrows(env, lp, 0, 1, inst->intersectionsLengths[i], &rhs2, "G", &izero, inst->intersections[i], inst->values, NULL, &name2);

    int branchingConstraint = CPXgetnumrows(env, lp)-1;


    if(verbose)
    {
        CPXwriteprob(env, lp, "Model3.lp", NULL);
    }

    if (CPXlpopt(env, lp)) // solve the problem
        print_error(" Problems on CPXlpopt2");

    CPXgetobjval(env, lp, &rightobj);
    
    delta2 = epsilon > rightobj-obj ? epsilon : rightobj-obj;

    CPXdelrows(env, lp, branchingConstraint, branchingConstraint);
    

    if(verbose)
    {
        CPXwriteprob(env, lp, "Model4.lp", NULL);
    }
    return delta1*delta2;
}



/**
 * Called inside the callback that performs strong branching. 
 * Cycles through all the branching constraints, for the ones that can be used it computes
 * the lp relaxation of the two childs.
 * It returns the index of the best constraint found and its product score.
 *
 * @param inst instance of the scp
 * @param rootSolution solution of the parent node
 * @param obj object value of parent node
 * @param cstat base of the parent node
 * @param rstat base of the parent node
 * @param ogNumCols original number of columns
 * @param threadNum index of the current thread
 * @param a contains the index of the best constraint found
 */

double computeConstraintsProductScore(instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, int ogNumCols, int threadNum, int* a)
{
    CPXENVptr env = inst->envs[threadNum];
    CPXLPptr lp = inst->lps[threadNum];

    double Max=0;
    double progressPercStep = 10.0/100;
    double progressPerc = progressPercStep;
    
    double currScore=0;
   
    // cycle on all the sets
    for(int i=0; i < inst->numIntersections; i++)
    {
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        for(int k = 0; k < inst->intersectionsLengths[i]; k++)
        {
            sum+=rootSolution[inst->intersections[i][k]];
        }

        if(sum>1e-05*inst->intersectionsLengths[i] && sum<1-1e-05*inst->intersectionsLengths[i]) // set covering
        {
            currScore = computeStrongBranchingOnConstraint(inst, obj, cstat, rstat, ogNumCols, threadNum, i);
        }

        if(Max<currScore)
        {
            // printf("new max for Constraints: %f, previous was: %f\n", currScore, Max);
            Max=currScore;
            *a=i;
        }
    }
    return Max;    
}



/**
 * Called inside the callback that performs strong branching. 
 * Cycles through all the branching constraints, for the ones that can be used it computes
 * the lp relaxation of the two childs.
 * It returns the index of the best constraint found and its product score.
 *
 * @param inst instance of the scp
 * @param rootSolution solution of the parent node
 * @param obj object value of parent node
 * @param cstat base of the parent node
 * @param rstat base of the parent node
 * @param threadNum index of the current thread
 * @param a contains the index of the best constraint found
 */

double computeConstraintsProductScoreContainingVar(instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, int ogNumCols, int bestVar, int threadNum, int* a)
{
    CPXENVptr env = inst->envs[threadNum];
    CPXLPptr lp = inst->lps[threadNum];

    double epsilon = epsilon;
    double Max=0;
    double progressPercStep = 10.0/100;
    double progressPerc = progressPercStep;
    
    double currScore=0;
   

   int *constraintCounter = inst->constraintCounter; 
    int **varConstrTable = inst->varConstrTable;
    int **intersections = inst->intersections;
    int *intersectionsLengths = inst->intersectionsLengths;
    
    int constCounter = constraintCounter[bestVar]; 
    //printf("Number of constraint containing the variable = %d\n",constCounter);
    // cycle the constraints that contain the best variable
    
    for (int i = 0; i < constCounter; i++)
    {
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        
        int constraintIndex = varConstrTable[bestVar][i];
        int *constraint = intersections[constraintIndex];
        int constraintLength = intersectionsLengths[constraintIndex];
        
        // printf("considering constraint :");
        // print_array_int(constraint, constraintLength);
        
        for (int k = 0; k < constraintLength; k++)
        {
            sum += rootSolution[constraint[k]];
        }
        
        if (sum > epsilon * constraintLength && sum < 1 - epsilon * constraintLength) // set covering
        {
            currScore = computeStrongBranchingOnConstraint(inst, obj, cstat, rstat, ogNumCols, threadNum, constraintIndex);
            //printf("currScore = %f\n", currScore);
            if (Max < currScore)
            {
                Max = currScore;
                *a = constraintIndex;
            }
        }       
    }
    return Max;    
}

/**
 * Called inside the callback that performs strong branching. 
 * Cycles through all the branching constraints, for the ones that can be used it computes
 * the lp relaxation of the two childs.
 * It returns the index of the best constraint found and its product score.
 *
 * @param env Cplex environment
 * @param lp Cplex problem 
 * @param inst instance of the scp
 * @param rootSolution solution of the parent node
 * @param obj object value of parent node
 * @param cstat base of the parent node
 * @param rstat base of the parent node
 * @param numCols number of variables i the parent node
 * @param indices array of indices used to restore bounds
 * @param lu specify the higher bound
 * @param originalbd specify the value of the higher bound
 * @param a contains the index of the best constraint found
 */ 
double computeVariableProductScore(instance* inst, double* rootSolution, double obj, int threadNum, int bdcnt, const int *nodebeg, const int* indices, const char *lu, const double *bd)
{
    
    CPXENVptr env = inst->envs[threadNum];
    CPXLPptr lp = inst->lps[threadNum];
    int status;
    //make a second copy of the lp for the branch1
    CPXLPptr branch1 = CPXcloneprob(env, lp, &status);
    //CPXwriteprob(env, lp, "Model1.lp", NULL);

    double epsilon = 0.001;
    double leftobj = 0;
    double rightobj = 0;

    double Max=0;
    int usebase = 1;

    int* cstat = inst->cstats[threadNum];
    int* rstat = inst->rstats[threadNum];
    double currScore=0;
    // cycle on all the sets
    // printf("Constraint %d usable, sum  = %f\n", i , sum);
    
    // print_array_int(indices, bdcnt);
    int cnt0 = nodebeg[1];
    int cnt1 = bdcnt-nodebeg[1];
   
    
    int indices0[cnt0];
    int indices1[cnt1];
    
    char lu0[cnt0];
    char lu1[cnt1];
    
    double bd0[cnt0];
    double bd1[cnt1];
    double bdReverse[cnt0];

    for(int i=0; i<bdcnt; i++)
    {
        if(i < cnt0)
        {
            indices0[i]=indices[i];
            lu0[i]=lu[i];
    
            bd0[i]=bd[i];

            if(bd[i]==0)
            {
                bdReverse[i]=1;
            }
            else
            {
                bdReverse[i]=0;
            }
        }
        else
        {
            indices1[i-cnt0]=indices[i];
            lu1[i-cnt0]=lu[i];
            bd1[i-cnt0]=bd[i];
        }
    }

    // test 0 branch
    if(usebase)
        CPXcopybase (env, lp, cstat, rstat);

    CPXchgbds(env, lp, cnt0, indices0, lu0, bd0);
 
    // CPXwriteprob(env, lp, "Model2.lp", NULL);

    if (CPXlpopt(env, lp)) // solve the problem
        print_error(" Problems on CPXlpopt1");

    CPXgetobjval(env, lp, &leftobj);
            
    double delta0 = epsilon > leftobj-obj ? epsilon : leftobj-obj;

    //reset problem
    // CPXchgbds(env, lp, cnt0, indices0, lu0, bdReverse);
    // CPXwriteprob(env, lp, "Model3.lp", NULL);

    //CPXwriteprob(env, branch1, "Model3.lp", NULL);



    // test 1 branch 
    // reset base
    if(usebase)
        CPXcopybase (env, branch1, cstat, rstat);
    
    CPXchgbds(env, branch1, cnt1, indices1, lu1, bd1);
    //CPXwriteprob(env, branch1, "Model4.lp", NULL);

    if (CPXlpopt(env, branch1)) // solve the problem
        print_error(" Problems on CPXlpopt2");

    CPXgetobjval(env, branch1, &rightobj);
    
    double delta1 = epsilon > rightobj-obj ? epsilon : rightobj-obj;

    currScore = delta0*delta1;
    //printf("variableProduct score = %f\n",currScore);
    return currScore;    
}


void saveNodeInfo(CPXCENVptr env, void *cbdata, int wherefrom, int threadNum, instance* inst)
{

    double epsilon = 0.0001;

    if(inst->nInfoIndex[threadNum] == inst->nInfoLengths[threadNum])
    {
        //printf("resizing info, currlen = %d, currindex = %d \n",inst->nInfoLengths[threadNum],inst->nInfoIndex[threadNum]);
        //resize nInfo
        inst->nInfoLengths[threadNum]*=2;
        inst->nInfo[threadNum] = (nodeInfo*) realloc(inst->nInfo[threadNum], inst->nInfoLengths[threadNum]*sizeof(nodeInfo));
    }

    int node;
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
    
    double value;
    CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &value);

    int depth;
    CPXgetcallbacknodeinfo(env, cbdata, wherefrom, 0, CPX_CALLBACK_INFO_NODE_DEPTH, &depth);
    
    double varFracPerc=0;
    double *x = (double *)malloc(inst->num_cols * sizeof(double));

    //printf("LP relaxation solved optimally it has the objective %f\n", obj);
    CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, inst->num_cols - 1);

    for(int i = 0; i<inst->num_cols; i++)
    {
        if (x[i]>epsilon && x[i]< 1 - epsilon)
        {
            varFracPerc++;
        }
    }
   
    varFracPerc = varFracPerc/inst->num_cols*100;
    int index = inst->nInfoIndex[threadNum];
    nodeInfo* nInfoThread = inst->nInfo[threadNum];
    
    nInfoThread[index].num = node;
    nInfoThread[index].value = value;
    nInfoThread[index].depth = depth;
    nInfoThread[index].varFracPerc = varFracPerc;
    printf("node = %d, value = %f, depth = %d, percentage of fractional variables = %f\n", inst->nInfo[threadNum][inst->nInfoIndex[threadNum]].num, inst->nInfo[threadNum][inst->nInfoIndex[threadNum]].value, inst->nInfo[threadNum][inst->nInfoIndex[threadNum]].depth, inst->nInfo[threadNum][inst->nInfoIndex[threadNum]].varFracPerc);

    inst->nInfoIndex[threadNum]++;
    free(x);
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
int legacyBranchingCallbackStrong(CPXCENVptr env, void *cbdata, int wherefrom, void *cbhandle, int type, int sos, int nodecnt, int bdcnt, const int *nodebeg,
                            const int *indices, const char *lu, const double *bd, const double *nodeest, int *useraction_p)
{

    int printCPLEXBranching=0;
    if(printCPLEXBranching)
    {
        printf("Number of branches = %d\n", nodecnt);
        if(nodecnt==2)
        {
            int cnt0 = nodebeg[1];
            int cnt1 = bdcnt-nodebeg[1];

            int indices0[cnt0];
            int indices1[cnt1];
            
            char lu0[cnt0];
            char lu1[cnt1];
            
            double bd0[cnt0];
            double bd1[cnt1];
            double bdReverse[cnt0];

            for(int i=0; i<bdcnt; i++)
            {
                if(i < cnt0)
                {
                    indices0[i]=indices[i];
                    lu0[i]=lu[i];
                    bd0[i]=bd[i];
                }
                else
                {
                    indices1[i-cnt0]=indices[i];
                    lu1[i-cnt0]=lu[i];
                    bd1[i-cnt0]=bd[i];
                }
            }
        printf("Branch 0:\n");
        print_array_int(indices0, cnt0);
        print_array_char(lu0, cnt0);
        print_array(bd0, cnt0);

        printf("\nBranch 1:\n");
        print_array_int(indices1, cnt1);
        print_array_char(lu1, cnt1);
        print_array(bd1, cnt1);

        // double obj;
        // CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);

        // int child1, child2;
        // CPXbranchcallbackbranchbds(env, cbdata, wherefrom, cnt0, indices0, lu0, bd0, obj, NULL, &child1);
        // CPXbranchcallbackbranchbds(env, cbdata, wherefrom, cnt1, indices1, lu1, bd1, obj, NULL, &child1);
        // return 0;
        }
    }   
    
    double start = second();
    double end;
    instance *inst = (instance *)cbhandle; // casting of cbhandle
    int threadNum;
    CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_MY_THREAD_NUM, &threadNum);

    // Default branching
    if (inst->branching == 2)
    {
        inst->defaultBranching[threadNum]++;
        end = second();
        inst->timeInCallback[threadNum] += end - start;
        saveNodeInfo(env, cbdata, wherefrom, threadNum, inst);
        return 0;
    }
    else
    {
        // printf("type = %d\n",type);
        // printf("SOS = %d\n",sos);
        if (nodecnt == 2)
        {
            int node;
            CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
            // sort based on the reduced costs
            int foundConstraint = 0;
            int i = -1;
            double obj;

            double *x=inst->xs[threadNum];
            
            CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);
            
            //printf("LP relaxation solved optimally it has the objective %f\n", obj);
            CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, inst->num_cols - 1);

            double objval;
            CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &objval);

            // printf("Objval relaxation at node = %f\n", objval);
            double Max = 0;
            
            // strong branching
            CPXLPptr nodelp;
            CPXgetcallbacknodelp(env, cbdata, wherefrom, &nodelp);

            int cur_numcols = CPXgetnumcols (env, nodelp);
            int cur_numrows = CPXgetnumrows (env, nodelp);
            getBase(env, nodelp, inst, threadNum, cur_numcols, cur_numrows);          
            
            int* cstat = inst->cstats[threadNum];
            int* rstat = inst->rstats[threadNum];
        
            int status;
            inst->lps[threadNum] = CPXcloneprob (inst->envs[threadNum], nodelp, &status);                
            // printf("clone status %d\n", status);

            // state that the problem is a lp problem
            CPXchgprobtype(inst->envs[threadNum], inst->lps[threadNum], 0);
            
            // char modelLog[50];
            // sprintf(modelLog, "ModelNode%d.lp", node);
            // printf("%s\n",modelLog);
            // CPXwriteprob(inst->envs[threadNum], nodeClone, modelLog, NULL);

            
            CPXgetub(inst->envs[threadNum], inst->lps[threadNum], inst->ogBd[threadNum], 0, cur_numcols-1);

            if(inst->constraintBranchVer == 1)
            {
                Max = computeConstraintsProductScoreContainingVar(inst, x, obj, cstat, rstat, cur_numcols, indices[0], threadNum, &i);
            }
            else
            {
                Max = computeConstraintsProductScore(inst, x, obj, cstat, rstat, cur_numcols, threadNum, &i);
            }
            
            double variableProductScore = computeVariableProductScore(inst, x, obj, threadNum, bdcnt, nodebeg, indices, lu, bd);

            // // get pseudocosts on the variables computed by CPLEX
            // double *pseudocostDown = inst->varPseudocostsDown[threadNum];
            // double *pseudocostUp = inst->varPseudocostsUp[threadNum];
            
            // CPXgetcallbackpseudocosts(env, cbdata, wherefrom, pseudocostUp, pseudocostDown, 0, inst->num_cols - 1);

            // double varCostDown = (0.001 > pseudocostDown[indices[0]] * x[indices[0]] ? 0.001 : pseudocostDown[indices[0]] * x[indices[0]]);
            // double varCostUp = (0.001 > pseudocostUp[indices[0]] * (1 - x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]] * (1 - x[indices[0]]));
            // double variableProductScore1 = varCostUp * varCostDown;
            // // printf("varproductScore = %f varproductscore1 = %f\n", variableProductScore, variableProductScore1);
            // // printf("comparing score for the constraint %f and score for the variable %f\n", Max, inst->delta * (varCostDown * varCostUp));
            
            
            if (i != -1 && Max > inst->delta * variableProductScore)
            {   
                // printf("Using branching constraint\n");
                // print_array_int(inst->intersections[i], inst->intersectionsLengths[i]);
                addBranchingChilds(env, cbdata, wherefrom, obj, i, threadNum, inst);
                *useraction_p = CPX_CALLBACK_SET;
                inst->constraintBranching[threadNum]++;
            }
            else
            {
                // printf("Using varaible branching on variable %d\n", indices[0]);
                inst->defaultBranching[threadNum]++;
            }
            end = second();
            inst->timeInCallback[threadNum] += end - start;
            saveNodeInfo(env, cbdata, wherefrom, threadNum, inst);
            return 0;
        }
        else
        {
            inst->defaultBranching[threadNum]++;
            end = second();
            inst->timeInCallback[threadNum] += end - start;
            saveNodeInfo(env, cbdata, wherefrom, threadNum, inst);
            return 0;
        }
    }
}


void sortIntersectionWRTreducedCosts(CPXCENVptr env, void *cbdata, int wherefrom, instance *inst)
{
    CPXLPptr nodelp;
    CPXgetcallbacknodelp (env, cbdata, wherefrom, &nodelp);
    double start = second();
  
    double* reducedCosts = (double *)calloc(inst->num_cols, sizeof(double));
   
    CPXgetdj(env, nodelp, reducedCosts, 0, inst->num_cols-1);
    
    // FILE *f;
    // f = fopen("ReducedCosts.txt", "w");
    // //fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
    // fprint_array(f, reducedCosts, inst->num_cols);
    // fclose(f);                    
    
    double end = second();

    printf("Retrieved reduced costs in %f up to now %f \n", end - start, end-inst->startTime);

    start = second();
    computeConstraintScoresReducedCosts(inst, reducedCosts);
    end = second();
    
    printf("Computed constraint score in %f up to now %f \n", end - start, end-inst->startTime);

        
    // f = fopen("ScoresConstr.txt", "w");
    // fprint_array(f, inst->constraintScoresD, inst->numIntersections);
    // fclose(f);        

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


    // f = fopen("ScoresConstrSorted1.txt", "w");
    // // fprint_array_int_int2(f, inst->intersections, inst->intersectionsLengths, inst->numIntersections);
    // // fprint_array_int(f, inst->intersectionsLengths, inst->numIntersections);
    // fprint_array(f, inst->constraintScoresD, inst->numIntersections);
    // fclose(f);        

    free(aux);
    free(aux1);
    free(aux2);
    free(reducedCosts);
}
          

void getBase(CPXCENVptr env, CPXLPptr lp, instance* inst, int threadNum, int cur_numcols, int cur_numrows)
{
    while(cur_numcols > inst->cstatDims[threadNum])
    {
        //resize cstat
        inst->cstatDims[threadNum]*=2;
        inst->cstats[threadNum] = (int*) realloc(inst->cstats[threadNum], inst->cstatDims[threadNum]*sizeof(int));

    }

    while(cur_numrows > inst->rstatDims[threadNum])
    {
        //resize rstat
        inst->rstatDims[threadNum]*=2;
        inst->rstats[threadNum] = (int*) realloc(inst->rstats[threadNum], inst->rstatDims[threadNum]*sizeof(int)); 
    }
    CPXgetbase (env, lp, inst->cstats[threadNum], inst->rstats[threadNum]);
}


/**
 * Save the details of the computations in a .csv file
 * @param inst instance of the scp
 * 
 */

void saveNodeInfoToFile(instance *inst)
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
    sprintf(folder, "../results/NodeInfo/%d_%d_%d_%.1f_%d_%d_%d_%d_%d", inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->reverse, inst->repeatedFirst, inst->sort, inst->average);
    sprintf(fileName, "../results/NodeInfo/%d_%d_%d_%.1f_%d_%d_%d_%d_%d/%s_%d.csv", inst->callback, inst->branching, inst->constraintBranchVer, inst->delta, inst->lookAhead, inst->reverse, inst->repeatedFirst, inst->sort, inst->average, input_file, inst->seed);
    
    
    printf("Saving results in file %s\n", fileName);
    // create folder if not present
    struct stat st = {0};
    if (stat(folder, &st) == -1) {
        mkdir(folder, 0700);
    }

    FILE *f;
    f = fopen(fileName, "w");
    fprintf(f, "Instance,Seed,NodeIndex,Gap,depth,percentageFracVars\n"); 
    
    for(int i=0; i<inst->threads; i++)
    {
        for(int j=0; j < inst->nInfoIndex[i]; j++)
        {
            fprintf(f, "%s,%d,%d,%f,%d,%f\n", inst->input_file, inst->seed, inst->nInfo[i][j].num, (abs(inst->nInfo[i][j].value - inst->bestVal))/inst->bestVal, inst->nInfo[i][j].depth, inst->nInfo[i][j].varFracPerc); 
        }    
    }
    fclose(f);
}
