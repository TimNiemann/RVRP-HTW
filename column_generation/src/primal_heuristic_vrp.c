/**@file   primal_heuristic_vrp.c
 * @brief  Functions for generating an inital solution
 * @author Lukas Schürmann, University Bonn
 */

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "scip/scip.h"

#include "tools_data.h"
#include "probdata_vrp.h"
#include "pricer_vrp.h"
#include "tools_vrp.h"
#include "labeling_algorithm_vrp.h"
#include "pricing_heuristic_vrp.h"


/** heuristic to calculate a good VRP solution based current LP-solution */
SCIP_RETCODE dispatchingHeuristic(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_PRICERDATA*    pricerdata
){
    model_data* modeldata;              // underlying modeldata
    tuple* valonday;                    // values of the day variables
    SCIP_VAR** modelvars;               // current variables
    SCIP_VARDATA* vardata;              // data of the variables
    SCIP_Real lpval;                    // value of the variable in the current solution
    int nvars;                          // number of variables in the current LP
    int i, j;
    int customer, day;
    int num_dayvars;
    double totalobj;
    SCIP_Bool isfeasible;
    SCIP_Bool localForbidden;
    char name[SCIP_MAXSTRLEN];
    char strtmp[SCIP_MAXSTRLEN];

    int duration;
    int nvisited;
    int* dayofcustomer;

    int** tour;                                 // current tour for each day
    int* tourlength;                            // current tour length for each day
    double* tourobj;                            // current tour objective value for each day

    solutionWindow*** solutionwindows;

    /* allocation & initialization */
    probdata = SCIPgetProbData(scip);
    modeldata = probdata->modeldata;
    num_dayvars = (modeldata->nC - 1) * modeldata->nDays;

    SCIP_CALL( SCIPallocMemoryArray( scip, &solutionwindows, modeldata->nDays));
    SCIP_CALL( SCIPallocMemoryArray( scip, &valonday, num_dayvars));
    SCIP_CALL( SCIPallocMemoryArray( scip, &dayofcustomer, modeldata->nC - 1));
    SCIP_CALL( SCIPallocMemoryArray( scip, &tour, modeldata->nDays));
    SCIP_CALL( SCIPallocMemoryArray( scip, &tourlength, modeldata->nDays));
    SCIP_CALL( SCIPallocMemoryArray( scip, &tourobj, modeldata->nDays));

    for(i = 0; i < modeldata->nDays; i++)
    {
        tourobj[i] = 0.0;
        tourlength[i] = 0;
        SCIP_CALL( SCIPallocMemoryArray( scip, &tour[i], modeldata->nC - 1));
    }

    /* we use a mapping from two-dimensional to one-dimensional vector */
    for(i = 0; i < modeldata->nC - 1; i++)
    {
        dayofcustomer[i] = -1;
        for(j = 0; j < modeldata->nDays; j++)
        {
            valonday[i + j * (modeldata->nC - 1)].value = 0.0;
            valonday[i + j * (modeldata->nC - 1)].index = i + j * (modeldata->nC - 1);
        }
    }
    SCIP_CALL( SCIPgetVarsData(scip, &modelvars, &nvars, NULL, NULL, NULL, NULL));

    /* calculate day-customer-variable values for current LP-solution */
    for(i = 0; i < nvars; i++)
    {
        lpval = SCIPvarGetLPSol(modelvars[i]);
        if(SCIPisSumPositive(scip, lpval))
        {
            vardata = SCIPvarGetData(modelvars[i]);
            if(vardata->tourlength == 0) continue;
            day = vardata->day;
            for(j = 0; j < vardata->tourlength; j++)
            {
                valonday[vardata->customertour[j] + day * (modeldata->nC - 1)].value += lpval;
            }
        }
    }
    /* sort after descending order */
    qsort(valonday, num_dayvars, sizeof(valonday[0]), cmp_vrp);
    nvisited = 0;
    /* try to add each customer to a fitting day based on day-variables */
    for(i = 0; i < num_dayvars; i++)
    {
        if(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) > 1 && SCIPisSumEQ(scip, valonday[i].value, 0.0))
        {
            break;
        }
        day = (int) (valonday[i].index / (modeldata->nC - 1));
        customer = valonday[i].index - day * (modeldata->nC - 1);
        if(dayofcustomer[customer] >= 0) continue;
        if(addNodeToTour(scip, modeldata, tour[day], &tourlength[day], &tourobj[day], customer, day, NULL, DBL_MAX))
        {
            dayofcustomer[customer] = day;
            nvisited++;
        }
        if(nvisited == modeldata->nC - 1) break;
    }
    if(nvisited == modeldata->nC - 1)
    {
        /* if every customer got assigned, improve the solution */
        SCIP_CALL( twoNodeShift(scip, modeldata, tour, tourlength, tourobj, dayofcustomer));
        totalobj = 0.0;
        for(day = 0; day < modeldata->nDays; day++)
        {
            totalobj += tourobj[day];
        }

        if(totalobj < SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)))
        {
            /* If this solution is better than the currently best one known, add new variables to SCIP */
            SCIP_SOL* solution;
            SCIP_CALL(SCIPcreateSol(scip, &solution, NULL));
            int newvars = 0;
            for(day = 0; day < modeldata->nDays; day++)
            {
                if (tourlength[day] > 0)
                {
                    newvars++;
                    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%2d: ", "primalDispatch", day);
                    for (i = 0; i < tourlength[day]; i++)
                    {
                        (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", tour[day][i]);
                        strcat(name, strtmp);
                    }
                    double obj = computeObjValue(scip, modeldata, &(solutionwindows[day]), &isfeasible, tour[day], &duration, tourlength[day], day);
                    assert(duration <= WORKTIME_LIMIT);
                    assert(obj == tourobj[day]);
                    assert(isfeasible);
                    SCIP_CALL( SCIPcreateColumn(scip, probdata, name, TRUE, tourobj[day], tour[day], tourlength[day], duration, solutionwindows[day], day));
                    SCIP_CALL( freeSolutionWindowArray(scip, solutionwindows[day], tourlength[day]) );
                    SCIP_CALL( SCIPsetSolVal(scip, solution, probdata->vars[probdata->nvars - 1], 1.0) );
                }else
                {
                    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "t_initEmpty_%2d", day);
                    assert(strcmp(name, SCIPvarGetName(probdata->vars[day])) == 0);
                    SCIP_CALL( SCIPsetSolVal(scip, solution, probdata->vars[day], 1.0) );
                }
            }
            SCIP_CALL( SCIPaddSol(scip, solution, &isfeasible));
            assert(isfeasible);
        }
    }

    for(day = 0; day < modeldata->nDays; day++)
    {
        SCIPfreeMemoryArray(scip, &tour[day]);
    }
    SCIPfreeMemoryArray(scip, &solutionwindows);
    SCIPfreeMemoryArray(scip, &tourobj);
    SCIPfreeMemoryArray(scip, &tourlength);
    SCIPfreeMemoryArray(scip, &tour);

    SCIPfreeMemoryArray(scip, &dayofcustomer);
    SCIPfreeMemoryArray(scip, &valonday);
    return SCIP_OKAY;
}