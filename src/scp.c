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
    if (inst->branching == 0)
    {
        inst->defaultBranching[threadNum]++;
        double end = second();
        inst->timeInCallback[threadNum] += end - start;
        return 0;
    }
    else
    {
        // only if CPLEX propose a branching that generates 2 new nodes.
        if (nodecnt == 2)
        {
            int foundConstraint = 0;
            int i = -1;
            double *x = (double *)calloc(inst->num_cols, sizeof(double));
            double obj;

            // int node;
            // CPXgetcallbackinfo(env, cbdata, wherefrom, CPX_CALLBACK_INFO_NODE_COUNT, &node);
            double *pseudocostDown = (double *)calloc(inst->num_cols, sizeof(double));
            double *pseudocostUp = (double *)calloc(inst->num_cols, sizeof(double));
            CPXgetcallbackpseudocosts(env, cbdata, wherefrom, pseudocostUp, pseudocostDown, 0, inst->num_cols - 1);

            CPXgetcallbacknodeobjval(env, cbdata, wherefrom, &obj);
            //printf("LP relaxation solved optimally it has the objective %f\n", obj);
            CPXgetcallbacknodex(env, cbdata, wherefrom, x, 0, inst->num_cols - 1);
            double Max = 0;
            double sec1 = second();
            int findBestMethod = 2; // 0 naive 1 smart 2 limited

            switch (findBestMethod)
            {
                case 0:
                {
                    Max = findBestBranchingConstraint(&i, x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
                case 1:
                {
                    Max = findBestBranchingConstraintSmart(&i, x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
                case 2:
                {
                    Max = findBestBranchingConstraintContainingVar(&i, indices[0], x, pseudocostDown, pseudocostUp, inst);
                    break;
                }
            }

            inst->timeFindingConstraint[threadNum] += second() - sec1;
            
            double varCostDown = (0.001 > pseudocostDown[indices[0]] * x[indices[0]] ? 0.001 : pseudocostDown[indices[0]] * x[indices[0]]);
            double varCostUp = (0.001 > pseudocostUp[indices[0]] * (1 - x[indices[0]]) ? 0.001 : pseudocostUp[indices[0]] * (1 - x[indices[0]]));

            if (i != -1 && Max > 2 * (varCostDown * varCostUp))
            {
                addBranchingChilds(env, cbdata, wherefrom, obj, i, inst);
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
 * @returns 0 if executed successfully
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

            //printf("Constraint %d  not usable\n",inst->interSetStart[i]+j);
            constraintScorePrevision[i] = epsilon * epsilon;
        }
        //print_array(productScoreConstraint, inst->interSetStart[i]+j+1);
        //printf("breakpoint\n");
        //printf("score: %f\n", productScoreConstraint[inst->interSetStart[i]+j]);
        if (Max < constraintScorePrevision[i])
        {
            Max = constraintScorePrevision[i];
            best_i = i;
            //printf("new max for Constraints prevision: %f\n", Max);
        }
    }
    *n = best_i;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return (Max);
}

double findBestBranchingConstraintContainingVar(int *n, int bestVar, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double *estimateScoreDown = (double *)calloc(inst->constraintCounter[bestVar], sizeof(double));
    double *estimateScoreUp = (double *)calloc(inst->constraintCounter[bestVar], sizeof(double));
    double *constraintScorePrevision = (double *)calloc(inst->constraintCounter[bestVar], sizeof(double));
    double epsilon = 0.0001;

    double Max = 0;
    int best_i = -1;

    // cycle the constraints that contain the best variable
    for (int i = 0; i < inst->constraintCounter[bestVar]; i++)
    {
        double pseudoDown;
        double pseudoUp;

        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        pseudoDown = 0;
        pseudoUp = INT_MAX;
        
        int constraintIndex = inst->varConstrTable[bestVar][i];
        int *constraint = inst->intersections[constraintIndex];    
        int constraintLength = inst->intersectionsLengths[constraintIndex];

        // printf("***Managing new constraint %d ***\n", inst->interSetStart[i]+j);
        for (int k = 0; k < constraintLength; k++)
        {
            // printf("variable %d has value %f and down pseudocost = %f\n", inst->intersections[inst->interSetStart[i]+j][k], rootSolution[inst->intersections[inst->interSetStart[i]+j][k]], pseudocostDown[inst->intersections[inst->interSetStart[i]+j][k]]);
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            int currVar = constraint[k];
            sum += x[currVar];
            pseudoDown += pseudocostDown[currVar] * x[currVar];
            // printf("current downpseudocost = %f\n", pseudoDown);
            double candidate = pseudocostUp[currVar];

            pseudoUp = (pseudoUp > candidate ? candidate : pseudoUp);
        }

        pseudoUp *= (1 - sum);

        if (sum > 1e-05 * constraintLength && sum < 1 - 1e-05 * constraintLength) // set covering
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

            //printf("Constraint %d  not usable\n",inst->interSetStart[i]+j);
            constraintScorePrevision[i] = epsilon * epsilon;
        }
        //print_array(productScoreConstraint, inst->interSetStart[i]+j+1);
        //printf("breakpoint\n");
        //printf("score: %f\n", productScoreConstraint[inst->interSetStart[i]+j]);
        if (Max < constraintScorePrevision[i])
        {
            Max = constraintScorePrevision[i];
            best_i = constraintIndex;
            //printf("new max for Constraints prevision: %f\n", Max);
        }
    }
    *n = best_i;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return (Max);
}



void solveUsingLegacyCallback(CPXENVptr env, CPXLPptr lp, instance *inst)
{
    char str[80];

    if (inst->branching == 0)
        sprintf(str, "logfile_defaultBranchinglegacy.txt");
    else
        sprintf(str, "logfile_constraintBranchinglegacy.txt");
    CPXsetlogfilename(env, str, "w");

    CPXsetbranchcallbackfunc(env, legacyBranchingCallback, inst);

    if (inst->threads == 0)
        CPXgetnumcores(env, &inst->threads);
    CPXsetintparam(env, CPX_PARAM_THREADS, inst->threads);

    CPXsetintparam(env, CPXPARAM_MIP_Strategy_CallbackReducedLP, CPX_OFF);
    // CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);//CPX_VARSEL_PSEUDO);//
    CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_PSEUDO);

    inst->constraintBranching = (int *)calloc(inst->threads, sizeof(int));
    inst->defaultBranching = (int *)calloc(inst->threads, sizeof(int));

    inst->timeInCallback = (double *)calloc(inst->threads, sizeof(double));

    inst->timeFindingConstraint = (double *)calloc(inst->threads, sizeof(double));

    if (inst->branching == 1)
    {
        //CPXsetintparam(env, CPXPARAM_MIP_Strategy_VariableSelect, CPX_VARSEL_STRONG);//CPX_VARSEL_PSEUDO);//

        // Get the rows of the problem
        int nnz;
        int *izero = (int *)calloc(inst->num_rows, sizeof(int));
        int *indexes = (int *)calloc(inst->num_rows * inst->num_cols, sizeof(int));         // this one is big!
        double *values = (double *)calloc(inst->num_rows * inst->num_cols, sizeof(double)); // this one is big too!
        int surplus_p;

        CPXgetrows(env, lp, &nnz, izero, indexes, values, inst->num_rows * inst->num_cols, &surplus_p, 0, inst->num_rows - 1);

        indexes = (int *)realloc(indexes, nnz * sizeof(int));     // realloc should free the unused spaced of the array
        values = (double *)realloc(values, nnz * sizeof(double)); // realloc should free the unused spaced of the array
        
        double start = second();
        populateIntersectionsOf2(izero, indexes, nnz, inst);
        double end = second();
        printf("Found %d constraint intersections sorted in %f\n", inst->numIntersections, end - start);

        // start = second();
        // populateIntersectionsOf2Sorted(izero, indexes, nnz, inst);
        // end = second();
        // printf("Found %d constraint intersections unsorted in %f\n", inst->numIntersections, end - start);

        // start = second();
        // populateIntersectionsOf2NoDup(izero, indexes, nnz, inst);
        // end = second();
        // printf("Found %d constraint intersections sorted without duplicates in %f\n", inst->numIntersections, end - start);


        // start = second();
        // computeVariableFrequency(izero, indexes, nnz, inst);
        // //populateIntersectionsOf4(inst);// not working yet
        // end = second();
        // printf("Computed variable frequencies in %f\n", end - start);

        start = second();
        populateVariableConstraintTable(inst);
        end = second();
        printf("Populated variable-constraint table in %f\n", end - start);

        free(izero);
        free(indexes);
        free(values);
    }

    if (CPXmipopt(env, lp)) // solve the problem
        print_error(" Problems on CPXmipopt");

    inst->solution = (double *)calloc(inst->num_cols, sizeof(double));
    double obj;

    if (CPXgetx(env, lp, inst->solution, 0, inst->num_cols - 1) == 0 && CPXgetbestobjval(env, lp, &obj) == 0)
        printf("Solution found with objval = %f\n", obj);
    else
        print_error("Problems in getting the solution");

    int constraintCount = 0;
    int defaultCount = 0;
    double callbackTime = 0;
    double findingConstraintTime = 0;

    if (inst->branching == 1)
    {
        for (int i = 0; i < inst->numIntersections; i++)
        {
            free(inst->intersections[i]);
        }
        free(inst->intersections);
        free(inst->variableFreq);
    }

    for (int i = 0; i < inst->threads; i++)
    {
        constraintCount += (inst->constraintBranching[i]);
        defaultCount += (inst->defaultBranching[i]);
        callbackTime += (inst->timeInCallback[i]);
        findingConstraintTime += (inst->timeFindingConstraint[i]);
    }

    saveComputationResults(inst);
    printf("%d constraint branching\n%d default branching\n", constraintCount, defaultCount);

    printf("Total time in callback = %f which is %f%% of the total\n", callbackTime, 100 * callbackTime / (second() - inst->startTime));
    printf("Time spent in callback choosing the best constraint: %f\n", findingConstraintTime);

    free(inst->solution);
    free(inst->constraintBranching);
    free(inst->defaultBranching);
}


double findBestBranchingConstraintSmart(int *n, double *x, double *pseudocostDown, double *pseudocostUp, instance *inst)
{
    double *estimateScoreDown = (double *)calloc(inst->numIntersections, sizeof(double));
    double *estimateScoreUp = (double *)calloc(inst->numIntersections, sizeof(double));
    // initialize estimateScoreUp to the max int value
    for (int i = 0; i < inst->numIntersections; i++)
    {
        estimateScoreUp[i] = INT_MAX;
    }

    double *constraintScorePrevision = (double *)calloc(inst->numIntersections, sizeof(double));

    double *sum = (double *)calloc(inst->numIntersections, sizeof(double));

    double epsilon = 0.0001;

    // cycle through all the variables
    for (int i = 0; i < inst->num_cols; i++)
    {
        if (x[i] > epsilon)
        {
            // cycle all the constraints in which the variable is present
            for (int j = 0; j < inst->constraintCounter[i]; j++)
            {
                sum[inst->varConstrTable[i][j]] += x[i];
                estimateScoreDown[inst->varConstrTable[i][j]] += pseudocostDown[i] * x[i];
                double candidate = pseudocostUp[i];
                estimateScoreUp[inst->varConstrTable[i][j]] = (estimateScoreUp[inst->varConstrTable[i][j]] > candidate ? candidate : estimateScoreUp[inst->varConstrTable[i][j]]);
            }
        }
        // else
        // {
        //    for(int j = 0; j < inst->constraintCounter[i]; j++)
        //    {
        //       double candidate = pseudocostUp[i];
        //       estimateScoreUp[inst->varConstrTable[i][j]] = (estimateScoreUp[inst->varConstrTable[i][j]] > candidate ? candidate : estimateScoreUp[inst->varConstrTable[i][j]]);
        //    }
        // }
    }

    double Max = 0;
    int best_i = -1;

    for (int i = 0; i < inst->numIntersections; i++)
    {
        // constraint is usable
        if (sum[i] > 1e-05 * (inst->intersectionsLengths[i]) && sum[i] < 1 - 1e-05 * (inst->intersectionsLengths[i]))
        {
            //printf("violated sum\n");
            // multiply the minimum pseudocode up for the fraction that is needed for the sum to reach 1
            estimateScoreUp[i] = estimateScoreUp[i] * (1 - sum[i]);
            constraintScorePrevision[i] = estimateScoreUp[i] * estimateScoreDown[i];
        }
        if (Max < constraintScorePrevision[i])
        {
            Max = constraintScorePrevision[i];
            best_i = i;
            //printf("new max for Constraints prevision: %f\n", Max);
        }
    }
    *n = best_i;
    free(estimateScoreDown);
    free(estimateScoreUp);
    free(constraintScorePrevision);
    return (Max);
}

void addBranchingChilds(CPXCENVptr env, void *cbdata, int wherefrom, double obj, int i, instance *inst)
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
    double *values = (double *)calloc(inst->intersectionsLengths[i], sizeof(double));
    for (int n = 0; n < inst->intersectionsLengths[i]; n++)
    {
        values[n] = 1;
    }

    /* We want to branch. */
    int child1, child2;
    // print_array_int(inst->intersections[inst->interSetStart[i]+j], inst->interSetLen[i]);
    // print_array(values, inst->interSetLen[i]);
    // printf("adding constraints with %d variables, rhs1=%f, rhs2=%f, sense1=%c, sense2=%c, izero=%d\n", inst->interSetLen[i], rhs1, rhs2, sense1, sense2, izero);

    if (CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->intersectionsLengths[i], &rhs1, "E", &izero, inst->intersections[i], values, obj, NULL, &child1) != 0)
    {
        print_error("First branch was not succesful\n");
    }

    if (CPXbranchcallbackbranchconstraints(env, cbdata, wherefrom, 1, inst->intersectionsLengths[i], &rhs2, "G", &izero, inst->intersections[i], values, obj, NULL, &child2) != 0)
    {
        print_error("Second branch was not succesful\n");
    }
    free(values);
}

/**
 * it saves other than the details of the computations:
 * execution time,
 * best int
 * best obj val
 * mip gap
 * nodes explored
 * nodes remaining
 * 
 */

void saveComputationResults(instance *inst)
{
;
}
    