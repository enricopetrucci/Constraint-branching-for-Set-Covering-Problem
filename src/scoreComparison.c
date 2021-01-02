#include "scp.h"
#include "constraintIntersection.h"
#include "scoreComparison.h"


// main function
int scoreComparison(instance *inst)
{
    int error; // error = 0 -> no error

    //instantiates the enviroment
    CPXENVptr env = CPXopenCPLEX(&error);
    //instantiates the problem 
    CPXLPptr lp = CPXcreateprob(env, &error, "scp");

    // set random seed
    CPXsetintparam(env, CPX_PARAM_RANDOMSEED, inst->seed);

    // set time limit based on value passed as command line argument
    CPXsetdblparam(env, CPX_PARAM_TILIM, inst->timelimit);
    
    //CPXsetlogfilename(env, "logScoreComparison.txt", "w");

    // import problem from file
    CPXreadcopyprob(env, lp, inst->input_file, "LP");

    if (inst->threads==0)
        CPXgetnumcores(env, &inst->threads);

    CPXsetintparam(env, CPX_PARAM_THREADS, inst->threads);


    inst->num_cols = CPXgetnumcols(env, lp);
    inst->num_rows = CPXgetnumrows(env, lp);
    printf("rows=%d, cols=%d\n", inst->num_rows, inst->num_cols);

    int varInd[inst->num_cols];
    char ctype[inst->num_cols];
    for(int i=0; i<inst->num_cols; i++)
    {
        varInd[i]=i;
        ctype[i]='C';
    }
    
    CPXchgctype(env, lp, inst->num_cols, varInd, ctype);

    // state that the problem is an lp problem
    CPXchgprobtype(env, lp, 0);
 
    CPXwriteprob(env, lp, "cModel.lp", NULL);

    //CPXwriteprob(env, lp, "model.lp", NULL);
    
    double* productScoreVariables = (double *) calloc(inst->num_cols, sizeof(double));

    double* pseudocostDown = (double *) calloc(inst->num_cols, sizeof(double));
    double* pseudocostUp = (double *) calloc(inst->num_cols, sizeof(double));
    
    double epsilon = 0.001;

        double* rootSolution = (double *) calloc(inst->num_cols, sizeof(double));
    double obj;

    if (CPXlpopt(env, lp)) // solve the problem
    print_error(" Problems on CPXlpopt");

    if(CPXgetx(env, lp, rootSolution, 0, inst->num_cols-1)==0 && CPXgetobjval(env, lp, &obj)==0)
        printf("Solution found with objval = %f\n", obj);
    else
        print_error("Problems in getting the solution");

    int cur_numcols = CPXgetnumcols (env, lp);
    int cur_numrows = CPXgetnumrows (env, lp);
    
    int* cstat = (int *) malloc (cur_numcols*sizeof(int));
    int* rstat = (int *) malloc (cur_numrows*sizeof(int));
    
    CPXgetbase (env, lp, cstat, rstat);


    computeVariablesProductScores(env, lp, inst, rootSolution, obj, cstat, rstat, productScoreVariables, pseudocostDown, pseudocostUp, epsilon);
    
    FILE* f;
    f = fopen("ScoresVariable.txt", "w");
    
    for(int i=0; i < inst->num_cols; i++)
    {
        fprintf(f, "var = %d, score = %f, pseudoDown = %f pseudoUp = %f\n", i, productScoreVariables[i], pseudocostDown[i], pseudocostUp[i]);
    }
    fclose(f);

    //print_array(productScoreVariables, inst->num_cols);
    
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
    populateIntersectionsOf2Sorted(izero, indexes, nnz, inst);
    double end = second();
    printf("Found %d constraint intersections in %f\n", inst->numIntersections, end-start);

    // start = second();
    // populateIntesectionsOf2NoDup(izero, indexes, nnz, inst);
    // end = second();
    // printf("Found %d constraint intersections in %f\n", inst->numIntersections, end-start);

    double* constraintScorePrevision = (double *) calloc(inst->numIntersections, sizeof(double));
    double* estimateScoreDown = (double *) calloc(inst->numIntersections, sizeof(double));
    double* estimateScoreUp = (double *) calloc(inst->numIntersections, sizeof(double));

    computePrevisionConstraintsScores(inst, rootSolution, pseudocostDown, pseudocostUp, constraintScorePrevision, estimateScoreDown, estimateScoreUp, epsilon, 0);

    f = fopen("ScoresEstimateConstraint.txt", "w");
    for(int i=0; i < inst->numIntersections; i++)
    {
        fprintf(f, "constraint = %d, score = %f, estimateScoreDown=%f, estimateScoreUp=%f\n", i, constraintScorePrevision[i], estimateScoreDown[i], estimateScoreUp[i]);
    }
    fclose(f);
    

    double* productScoreConstraints = (double *) calloc(inst->numIntersections, sizeof(double));
    double* scoreDown = (double *) calloc(inst->numIntersections, sizeof(double));
    double* scoreUp = (double *) calloc(inst->numIntersections, sizeof(double));


    computeConstraintsProductScores(env, lp, inst, rootSolution, obj, cstat, rstat, productScoreConstraints, scoreDown, scoreUp, epsilon);

    f = fopen("ScoresConstraint.txt", "w");
    for(int i=0; i < inst->numIntersections; i++)
    {
        fprintf(f, "constraint = %d, score = %f, scoreDown=%f, scoreUp=%f\n", i, productScoreConstraints[i], scoreDown[i], scoreUp[i]);
    }
    fclose(f);


    double MaxVarScore = 0;
    for(int i = 0; i < inst->num_cols; i++)
    {
        if(MaxVarScore<productScoreVariables[i])
        {
            MaxVarScore=productScoreVariables[i];
        }
    }

    double MaxConsScore = 0;
    for(int i = 0; i < inst->numIntersections; i++)
    {
        if(MaxConsScore<productScoreConstraints[i])
        {
            MaxConsScore=productScoreConstraints[i];
        }
    }
    f = fopen("ConstraintScores.txt", "w");
    fprint_array(f, productScoreConstraints, inst->numIntersections);
    fclose(f);
    
    f = fopen("VariableScores.txt", "w");
    fprint_array(f, productScoreVariables, inst->num_cols);
    fclose(f);
    
    printf("MaxVarScore = %f, MaxConsScore = %f\n", MaxVarScore, MaxConsScore);

    double MSEScore = 0;
    double MSEScoreDown = 0;
    double MSEScoreUp = 0;
    
    for(int i = 0; i<inst->numIntersections; i++)
    {
        MSEScore += pow(productScoreConstraints[i]-constraintScorePrevision[i], 2);
        MSEScoreDown += pow(scoreDown[i]-estimateScoreDown[i], 2);
        MSEScoreUp += pow(scoreUp[i]-estimateScoreUp[i], 2);
    }
    printf("Policy: Min. MSEScore = %f, MSEScoreDown = %f, MSEScoreUp = %f\n", MSEScore/inst->numIntersections, MSEScoreDown/inst->numIntersections, MSEScoreUp/inst->numIntersections);

    computePrevisionConstraintsScores(inst, rootSolution, pseudocostDown, pseudocostUp, constraintScorePrevision, estimateScoreDown, estimateScoreUp, epsilon, 1);

    MSEScore = 0;
    MSEScoreDown = 0;
    MSEScoreUp = 0;
    
    for(int i = 0; i<inst->numIntersections; i++)
    {
        MSEScore += pow(productScoreConstraints[i]-constraintScorePrevision[i], 2);
        MSEScoreDown += pow(scoreDown[i]-estimateScoreDown[i], 2);
        MSEScoreUp += pow(scoreUp[i]-estimateScoreUp[i], 2);
    }
    printf("Policy: Max. MSEScore = %f, MSEScoreDown = %f, MSEScoreUp = %f\n", MSEScore/inst->numIntersections, MSEScoreDown/inst->numIntersections, MSEScoreUp/inst->numIntersections);

    computePrevisionConstraintsScores(inst, rootSolution, pseudocostDown, pseudocostUp, constraintScorePrevision, estimateScoreDown, estimateScoreUp, epsilon, 2);

    MSEScore = 0;
    MSEScoreDown = 0;
    MSEScoreUp = 0;
    
    for(int i = 0; i<inst->numIntersections; i++)
    {
        MSEScore += pow(productScoreConstraints[i]-constraintScorePrevision[i], 2);
        MSEScoreDown += pow(scoreDown[i]-estimateScoreDown[i], 2);
        MSEScoreUp += pow(scoreUp[i]-estimateScoreUp[i], 2);
    }
    printf("Policy: Max. MSEScore = %f, MSEScoreDown = %f, MSEScoreUp = %f\n", MSEScore/inst->numIntersections, MSEScoreDown/inst->numIntersections, MSEScoreUp/inst->numIntersections);


    // free and close CPLEX model
    CPXfreeprob(env, &lp);
    CPXcloseCPLEX(&env); 
}



void computeVariablesProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, double* productScoreVariables, double* pseudocostDown, double* pseudocostUp, double epsilon)
{
    double* leftSolution = (double *) calloc(inst->num_cols, sizeof(double));
    double leftobj;
    
    double* rightSolution = (double *) calloc(inst->num_cols, sizeof(double));
    double rightobj;

    
    FILE *f;

    f = fopen("Solutions.txt", "w");

    fprint_array(f, rootSolution, inst->num_cols);
    double Max=0;

    for(int i=0; i<inst->num_cols; i++)
    {
        double zero = 0;
        double one = 1;
        double delta1;
        double delta2;
        if(rootSolution[i]!=0)
        {
            // first node
            CPXchgbds( env, lp, 1, &i, "B", &zero);
            CPXcopybase (env, lp, cstat, rstat);    
            if (CPXlpopt(env, lp)) // solve the problem
                print_error(" Problems on CPXmipopt");

           
            CPXgetx(env, lp, leftSolution, 0, inst->num_cols-1);
            CPXgetobjval(env, lp, &leftobj);
            //printf("Left Solution found with objval = %f\n", leftobj);
            
            fprint_array(f, leftSolution, inst->num_cols);
            delta1 = leftobj-obj;
            pseudocostDown[i]=delta1/rootSolution[i];
        }
        else
        {
            fprint_array(f, rootSolution, inst->num_cols);
            delta1 = epsilon;
        }

        if(rootSolution[i]!=1)
        {
            // second node        
            CPXchgbds( env, lp, 1, &i, "B", &one);
            CPXcopybase (env, lp, cstat, rstat);
            if (CPXlpopt(env, lp)) // solve the problem
                print_error(" Problems on CPXmipopt");

            CPXgetx(env, lp, rightSolution, 0, inst->num_cols-1);
            CPXgetobjval(env, lp, &rightobj);
            //printf("Right Solution found with objval = %f\n", rightobj);
            
            fprint_array(f, rightSolution, inst->num_cols);
            delta2 = rightobj-obj;

            pseudocostUp[i]=delta2/(1-rootSolution[i]);
        }
        else
        {
            fprint_array(f, rootSolution, inst->num_cols);
            delta2 = epsilon;
        }

        
        productScoreVariables[i] = delta1*delta2;
        
        // printf("Variable = %d, Score = %f, Delta1 = %f, Delta2= %f\n", i, productScoreVariables[i], delta1, delta2);
        CPXchgbds( env, lp, 1, &i, "L", &zero);
        CPXchgbds( env, lp, 1, &i, "U", &one);
        if(Max<productScoreVariables[i])
        {
            Max=productScoreVariables[i];
            printf("new max for variables: %f\n", Max);
        }
    }
    fclose(f);  
}


void computeConstraintsProductScores(CPXENVptr env, CPXLPptr lp, instance* inst, double* rootSolution, double obj, int* cstat, int* rstat, double* productScoreConstraints, double* scoreDown, double*scoreUp, double epsilon)
{
    double* leftSolution = (double *) calloc(inst->num_cols, sizeof(double));
    double leftobj;
    
    double* rightSolution = (double *) calloc(inst->num_cols, sizeof(double));
    double rightobj;

    //FILE *f;

    //f = fopen("SolutionsConstraints.txt", "w");

    //fprint_array(f, rootSolution, inst->num_cols);
    double Max=0;
    int branchingConstraint;
    double progressPercStep = 10.0/100;
    double progressPerc = progressPercStep;
    
    // cycle on all the sets
    for(int i=0; i < inst->numIntersections; i++)
    {
        double delta1;
        double delta2;
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        for(int k = 0; k < inst->intersectionsLengths[i]; k++)
        {
            sum+=rootSolution[inst->intersections[i][k]];
        }

        if(sum>1e-05*inst->intersectionsLengths[i] && sum<1-1e-05*inst->intersectionsLengths[i]) // set covering
        {
            //printf("Constraint %d usable, sum  = %f\n",inst->interSetStart[i]+j, sum);
            double rhs1 = 0;
            double rhs2 = 1;
            char sense1 = 'E';
            char sense2 = 'G';
            char* name1 = "down";
            char* name2 = "up";

            int izero = 0;
            // generating the array of values, all ones for the variables selected
            double* values = (double *) calloc(inst->intersectionsLengths[i], sizeof(double));
            for(int n=0; n<inst->intersectionsLengths[i]; n++)
            {
                values[n]=1;
            }
            
            CPXaddrows(env, lp, 0, 1, inst->intersectionsLengths[i], &rhs1, "E", &izero, inst->intersections[i], values, NULL, &name1);
            
            branchingConstraint = CPXgetnumrows(env, lp)-1;
            CPXcopybase (env, lp, cstat, rstat);
            if (CPXlpopt(env, lp)) // solve the problem
                print_error(" Problems on CPXlpopt1");

            CPXgetx(env, lp, leftSolution, 0, inst->num_cols-1);
            CPXgetobjval(env, lp, &leftobj);
            //printf("Left Solution found with objval = %f\n", leftobj);
            
            // fprint_array(f, leftSolution, inst->num_cols);
    
            delta1 = epsilon > leftobj-obj ? epsilon : leftobj-obj;
    
            CPXdelrows(env, lp, branchingConstraint, branchingConstraint);

            CPXaddrows(env, lp, 0, 1, inst->intersectionsLengths[i], &rhs2, "G", &izero, inst->intersections[i], values, NULL, &name2);

            branchingConstraint = CPXgetnumrows(env, lp)-1;
            CPXcopybase (env, lp, cstat, rstat);
            if (CPXlpopt(env, lp)) // solve the problem
                print_error(" Problems on CPXlpopt2");

            CPXgetx(env, lp, rightSolution, 0, inst->num_cols-1);
            CPXgetobjval(env, lp, &rightobj);
            //printf("Right Solution found with objval = %f\n", rightobj);
            
            // fprint_array(f, rightSolution, inst->num_cols);

            delta2 = epsilon > rightobj-obj ? epsilon : rightobj-obj;
    
            CPXdelrows(env, lp, branchingConstraint, branchingConstraint);
            
            free(values);
            scoreDown[i] = delta1;
            scoreUp[i] = delta2;
            productScoreConstraints[i] = delta1*delta2;
        }
        else
        {
            //printf("Constraint %d  not usable\n",inst->interSetStart[i]+j);
            productScoreConstraints[i] = epsilon*epsilon;
        }
        //print_array(productScoreConstraint, inst->interSetStart[i]+j+1);
        //printf("breakpoint\n");
        //printf("score: %f\n", productScoreConstraint[inst->interSetStart[i]+j]);
        if(Max<productScoreConstraints[i])
        {
            Max=productScoreConstraints[i];
            //printf("new max for Constraints: %f\n", Max);
        }
        //printf("counter = %d, next progress step = %d, or %f. ProgressPerc = %f, NumIntersection = %d \n", inst->interSetStart[i]+j, (int)(inst->numIntersections*progressPerc), (inst->numIntersections*progressPerc), progressPerc, inst->numIntersections);
        if(i == (int)(progressPerc * inst->numIntersections))
        {
            printf("Completed %.0f%% \n", progressPerc*100);
            progressPerc+=progressPercStep;
        }
    }
    //fclose(f);
}




//policy: 0 = min 1 = max 2 = mean
    
// try considering the number of zeros. compute MSE
void computePrevisionConstraintsScores(instance* inst, double* rootSolution, double* pseudocostDown, double* pseudocostUp, double* constraintScorePrevision, double* estimateScoreDown, double* estimateScoreUp, double epsilon, int policy)
{
    double Max = 0;
    for(int i=0; i < inst->numIntersections; i++)
    {
        double pseudoDown;
        double pseudoUp;
        // printf("checking intersection of lenght = %d\n", inst->interSetLen[i]);
        // printf("j from 0 to %d\n", (i < inst->numInterSet-1) ? (inst->interSetStart[i+1] - inst->interSetStart[i]) : inst->numIntersections - inst->interSetStart[i]);
        // cycle on all the intersection having the same lenght
        double sum = 0;
        // cycle all the variables in the intersection and get the sum
        pseudoDown=0;
        if(policy==0)
        {
            pseudoUp=INT_MAX; 
        }
        else
        {
            pseudoUp=0;
        }
        // printf("***Managing new constraint %d ***\n", inst->interSetStart[i]+j);
        for(int k = 0; k < inst->intersectionsLengths[i]; k++)
        {
            // printf("variable %d has value %f and down pseudocost = %f\n", inst->intersections[inst->interSetStart[i]+j][k], rootSolution[inst->intersections[inst->interSetStart[i]+j][k]], pseudocostDown[inst->intersections[inst->interSetStart[i]+j][k]]);
            //printf("k=%d index = %d\n", k, inst->intersections[inst->interSetStart[i]+j][k]);
            // printf("number %d, variable %d has coefficient %.20f\n",k, inst->intersections[inst->interSetStart[i]+j][k],x[inst->intersections[inst->interSetStart[i]+j][k]]);
            sum+=rootSolution[inst->intersections[i][k]];
            pseudoDown += pseudocostDown[inst->intersections[i][k]]*rootSolution[inst->intersections[i][k]];
            // printf("current downpseudocost = %f\n", pseudoDown);
            double candidate = pseudocostUp[inst->intersections[i][k]];
            
            if(policy == 0)
            {
                pseudoUp = (pseudoUp > candidate ? candidate : pseudoUp);
            }
            else if(policy ==  1)
            {
                pseudoUp = (pseudoUp < candidate ? candidate : pseudoUp);
            }
            else if(policy == 2)
            {
                pseudoUp += candidate;                    
            }
        }
        if(policy==2)
        {
            pseudoUp /= inst->intersectionsLengths[i];
        }

        pseudoUp *= (1-sum);

        if(sum>1e-05*inst->intersectionsLengths[i] && sum<1-1e-05*inst->intersectionsLengths[i]) // set covering
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
            //printf("new max for Constraints prevision: %f\n", Max);
        }
    }
}
    