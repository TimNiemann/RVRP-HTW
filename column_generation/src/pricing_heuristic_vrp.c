/**@file   pricing_heurisitc_vrp.c
 * @brief  heuristic pricing as an alternative to the labeling algorithmn
 * @author Tim Niemann, TU Braunschweig
 * @author Lukas Schürmann, University Bonn
 */

#include <assert.h>
#include <string.h>
#include <stdlib.h>

#include "scip/scip.h"
#include "tools_data.h"
#include "pricer_vrp.h"
#include "probdata_vrp.h"
#include "vardata_vrp.h"
#include "tools_vrp.h"
#include "scip/cons_setppc.h"
#include "pricing_heuristic_vrp.h"

/** tries to add customers to tour in order to lower the reduced cost */
static
SCIP_RETCODE extendColumn(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_Bool           isFarkas,
        int*                tour,
        int*                length,
        double*             obj,
        double*             sumdual,
        double*             dualvalues,
        SCIP_Bool*          isinTour,
        tuple*              tocheck,
        int                 ntocheck,
        int                 day
){
    int i, k;
    double threshold;

    for(k = 0; k < ntocheck; k++)
    {
        i = tocheck[k].index;
        if(!isinTour[i])
        {
            if(isFarkas && !SCIPisSumPositive(scip, dualvalues[i]))
                continue;
            threshold = isFarkas ? SCIP_DEFAULT_INFINITY : dualvalues[i];
            if(addNodeToTour(scip, probdata->modeldata, tour, length, obj, i, day, NULL, threshold))
            {
                *sumdual += dualvalues[i];
                isinTour[i] = TRUE;
            }
        }
    }
    return SCIP_OKAY;
}

/** tries to delete customers from a tour in order to lower the reduced cost */
static
SCIP_RETCODE decreaseColumn(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_Bool           isFarkas,
        int*                tour,
        int*                length,
        double*             obj,
        double*             sumdual,
        double*             dualvalues,
        SCIP_Bool*          isinTour,
        tuple*              tocheck,
        int                 ntocheck,
        int                 day
){
    int i, j, k, l, pos;
    double tmpobj;
    int tmplength;
    SCIP_Bool isfeasible;
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    assert(pricerdata != NULL);
    assert(pricerdata->isForbidden != NULL);

    for(l = ntocheck - 1; l >= 0; l--)
    {
        i = tocheck[l].index;
        if(*length <= 1) break;
        if(isinTour[i])
        {
            /* try to delete every customer that gets visited by the tour */
            if(isFarkas)
            {
                if(SCIPisSumNegative(scip, dualvalues[i]))
                {
                    k = 0;
                    /* delete customer from tour */
                    for(j = 0; j < *length; j++)
                    {
                        if(k == 1)
                        {
                            tour[j - 1] = tour[j];
                        }
                        else if(tour[j] == i)
                        {
                            if(isDeleteAllowed(tour, *length, probdata->modeldata->nC - 1, j, pricerdata->isForbidden))
                            {
                                pos = -1;
                                break;
                            }
                            k = 1;
                            pos = j;
                            continue;
                        }
                    }
                    if(pos == -1) continue;
                    /* update values */
                    *length -= 1;
                    *sumdual -= dualvalues[i];
                    isinTour[i] = FALSE;
                    *obj = computeObjValue(scip, probdata->modeldata, NULL, &isfeasible, tour, NULL, *length, day);

                    assert(isfeasible);
                }
            }else{
                k = 0;
                /* create tour without customer i */
                for(j = 0; j < *length; j++)
                {
                    if(k == 1)
                    {
                        tour[j - 1] = tour[j];
                    }
                    else if(tour[j] == i)
                    {
                        if(!isDeleteAllowed(tour, *length, probdata->modeldata->nC - 1, j, pricerdata->isForbidden))
                        {
                            pos = -1;
                            break;
                        }
                        k = 1;
                        pos = j;
                        continue;
                    }
                }
                if(pos == -1) continue;
                tmplength = *length - 1;
                /* calculate new objective value and check if feasible */
                tmpobj = computeObjValue(scip, probdata->modeldata, NULL, &isfeasible, tour, NULL, tmplength, day);
                if(!isfeasible){
                    /* undo changes */
                    for(j = *length - 1; j > pos; j--)
                    {
                        tour[j] = tour[j - 1];
                    }
                    tour[pos] = i;
                    continue;
                }
                /* update tour data if improvement has been found */
                if(SCIPisSumNegative(scip, - (*obj - tmpobj - dualvalues[i])))
                {
                    *length = tmplength;
                    *obj = tmpobj;
                    *sumdual -= dualvalues[i];
                    isinTour[i] = FALSE;
                }else{
                    /* undo changes */
                    for(j = *length - 1; j > pos; j--)
                    {
                        tour[j] = tour[j - 1];
                    }
                    tour[pos] = i;
                }
            }
        }
    }
    return SCIP_OKAY;
}

/** tries to exchange customers in tour by others in order to lower the reduced cost */
static
SCIP_RETCODE shiftColumn(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_Bool           isFarkas,
        int*                tour,
        int                 length,
        double*             obj,
        double*             sumdual,
        double*             dualvalues,
        SCIP_Bool*          isinTour,
        tuple*              tocheck,
        int                 ntocheck,
        int                 day
){
    int node, currentnode, i;
    int pos;
    double tmpobj;
    SCIP_Bool isfeasible;
    model_data* modeldata;
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    assert(pricerdata != NULL);
    assert(pricerdata->isForbidden != NULL);

    modeldata = probdata->modeldata;

    assert(length > 0);
    assert(day >= 0 && modeldata->nDays > day);

    for(i = 0; i < ntocheck; i++)
    {
        currentnode = tocheck[i].index;
        assert(currentnode >= 0 && currentnode < modeldata->nC - 1);
        if(!isinTour[currentnode])
        {
            /* try each unvisited customer on each spot in tour */
            for (pos = 0; pos < length; pos++)
            {
                if(!isExchangeAllowed(tour, length, modeldata->nC - 1, pos, currentnode, pricerdata->isForbidden))
                {
                    continue;
                }
                node = tour[pos];
                assert(0 <= node && node < modeldata->nC - 1);
                tour[pos] = currentnode;
                tmpobj = computeObjValue(scip, probdata->modeldata, NULL, &isfeasible, tour, NULL, length, day);
                if (!isfeasible)
                {
                    tour[pos] = node;
                    continue;
                }
                /* update tour data if improvement was found */
                if(SCIPisSumNegative(scip, - ((!isFarkas) * (*obj - tmpobj) + dualvalues[currentnode] - dualvalues[node])))
                {
                    *obj = tmpobj;
                    *sumdual += dualvalues[currentnode] - dualvalues[node];
                    isinTour[currentnode] = TRUE;
                    isinTour[node] = FALSE;
                }else{
                    /* else undo change */
                    tour[pos] = node;
                }
                break;
            }
        }
    }
    return SCIP_OKAY;
}

/** check for local improvements for a given variable */
static
SCIP_Bool investigateColumn(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_Bool           isFarkas,
        SCIP_Bool*          solfoundforday,
        SCIP_VAR*           var,
        double*             dualvalues,
        int*                newtour,
        int*                newtourlength,
        double*             newtourobj,
        char*               allnames
){
    int i, j;
    int day;
    model_data* modeldata;
    SCIP_VARDATA* vardata = NULL;
    SCIP_Bool* isinTour;                        // [node] = True iff node is contained in the current tour
    SCIP_Bool improvement;
    char name[SCIP_MAXSTRLEN];                  // name of the constraint
    char strtmp[SCIP_MAXSTRLEN];

    int* tour;                                  // current tour
    int tourLength;                             // current tour length
    double obj;                                 // current tour objective value
    double lhs;                                 // value of the left hand side of the corresponding dual constraint

    tuple* tocheck;                             // Active nodes of the day - sorted by value of the dual variables
    int ntocheck;

    modeldata = probdata->modeldata;

    improvement = FALSE;
    /* get var data */
    vardata = SCIPvarGetData(var);
    day = vardata->day;
    assert(day >= 0 && day < modeldata->nDays);
    if(solfoundforday[day]) // only max. one new column per day
    {
        return FALSE;
    }
    obj = SCIPvarGetObj(var);
    lhs = dualvalues[modeldata->nC - 1 + day];
    tourLength = vardata->tourlength;

    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tour, modeldata->nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &isinTour, modeldata->nC - 1) );

    for(i = 0; i < modeldata->nC - 1; i++)
    {
        isinTour[i] = FALSE;
    }
    for (j = 0; j < tourLength; j++)
    {
        tour[j] = vardata->customertour[j];
        lhs += dualvalues[vardata->customertour[j]];
        isinTour[vardata->customertour[j]] = TRUE;
        assert(tour[j] != modeldata->nC - 1);
    }

    /* sort active nodes of the day */
    ntocheck = getNumberOfNeighbors(scip, modeldata, modeldata->nC - 1, day);
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tocheck, ntocheck) );
    SCIP_CALL( sortNeighborsOfNode(scip, modeldata, day, modeldata->nC - 1, tocheck, ntocheck, dualvalues));
    /*********** local search: ***********/

    /* (i) try to add nodes to the tour */
    SCIP_CALL(extendColumn(scip, probdata, isFarkas, tour, &tourLength, &obj, &lhs, dualvalues, isinTour, tocheck,
                           ntocheck, day));
    if (tourLength > 0) {
        /* (ii) try to exchange nodes of the tour by unused ones */
        SCIP_CALL(shiftColumn(scip, probdata, isFarkas, tour, tourLength, &obj, &lhs, dualvalues, isinTour, tocheck,
                              ntocheck, day));
        /* (iii) try to delete nodes from the tour (seems counterproductive when searching for a feasible solution) */
        if (!isFarkas) {
            SCIP_CALL(decreaseColumn(scip, probdata, isFarkas, tour, &tourLength, &obj, &lhs, dualvalues, isinTour,
                                     tocheck, ntocheck, day));
        }
    }
    SCIP_CALL(rearrangeTour(scip, modeldata, tour, tourLength, &obj, day));

    SCIPfreeBlockMemoryArray(scip, &tocheck, ntocheck);

    /* check if there has been an improvement, i.e. if the manipulated tour has negative reduced costs */
    if(SCIPisSumPositive(scip, lhs - !isFarkas * obj)) // check for violation of dual constraint: lhs <= 0.0
    {
        improvement = TRUE;
        if(tourLength == 0){
            improvement = FALSE;
        }
    }

    if(improvement) // tour with negative reduced costs has been found!
    {
        assert(tourLength > 0);
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingLocSearch_%2d: ", day);
        for (i = 0; i < tourLength; i++)
        {
            (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", tour[i]);
            strcat(name, strtmp);
        }
        if (SCIPprobdataContainsVar(probdata, name)) // check if corresponding column already exists - happens due to aging
        {
            return FALSE;
        }
        /* safe new tour */
        *newtourobj = obj;
        *newtourlength = tourLength;
        for(i = 0; i < tourLength; i++)
        {
            newtour[i] = tour[i];
        }
        newtour[*newtourlength] = day;
        (void) SCIPsnprintf(allnames, SCIP_MAXSTRLEN, "%s", name);
        /* set day to used */
        solfoundforday[day] = TRUE;
    }

    SCIPfreeBlockMemoryArray(scip, &isinTour, modeldata->nC - 1);
    SCIPfreeBlockMemoryArray(scip, &tour, modeldata->nC - 1);
    return improvement;
}

/** local search pricing to find negative reduced cost tours */
static
SCIP_RETCODE localSearchPricing(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_Bool           isFarkas
){
    /* auxiliary variables */
    int i, j, k;
    int nvars, nbinvars;
    int addedCols;
    int duration;
    double tmpobj;
    SCIP_Bool isfeasible;

    int nsortedvars;
    tuple* sortedvars;              // sorted list of the linear primal variables != 0
    tuple* backupvars;

    int** newtours;                 // created tours with negative reduced costs
    int* newtourslength;
    double* newtoursobj;
    char allnames[MAX_HEURISTIC_TOURS][SCIP_MAXSTRLEN];

    double* dualvalues;             // current values of the dual variables
    SCIP_Bool* solfoundforday;      // We only want to find at most one new column per day

    SCIP_VAR** vars = NULL;
    solutionWindow** solutionwindows;
    model_data* modeldata;

    modeldata = probdata->modeldata;

    /* allocation & initialization */
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &dualvalues, modeldata->nC - 1 + modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &solfoundforday, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newtours, MAX_HEURISTIC_TOURS) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newtourslength, MAX_HEURISTIC_TOURS) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newtoursobj, MAX_HEURISTIC_TOURS) );
    for(i = 0; i < MAX_HEURISTIC_TOURS; i++)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &newtours[i], modeldata->nC) );
    }
    for(i = 0; i < modeldata->nDays; i++)
    {
        solfoundforday[i] = FALSE;
    }

    SCIP_CALL( getDualValues(scip, dualvalues, isFarkas));

    /* get current variables of the RMP */
    SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, &nbinvars, NULL, NULL, NULL));
    assert(nvars == nbinvars);
    /* sort these variables by their LP-value. Only consider variables with LP-value > 0 */
    nsortedvars = 0;
    for(i = 0; i < nvars; i++)
    {
        if(SCIPisSumPositive(scip, SCIPvarGetLPSol(vars[i])))
        {
            nsortedvars++;
        }
    }
    j = 0, k = 0;
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sortedvars, nsortedvars));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &backupvars, nvars - nsortedvars));
    for(i = 0; i < nvars; i++)
    {
        if(SCIPisSumPositive(scip, SCIPvarGetLPSol(vars[i])))
        {
            sortedvars[j].index = i;
            sortedvars[j++].value = SCIPvarGetLPSol(vars[i]);
        }else{
            backupvars[k].index = i;
            backupvars[k++].value = -SCIPgetVarRedcost(scip, vars[i]);
        }
    }
    qsort(sortedvars, nsortedvars, sizeof(sortedvars[0]), cmp_vrp);
    qsort(backupvars, nvars - nsortedvars, sizeof(backupvars[0]), cmp_vrp);

    addedCols = 0;          // number of found tours with negative reduced costs
    for(j = 0; j < nsortedvars; j++)    // start from the variable with the highest LP-value
    {
        if(addedCols == MAX_HEURISTIC_TOURS) break;     // stop if enough new tours have been found

        assert(vars[sortedvars[j].index] != NULL);

        /** start local search for variable */
        if(investigateColumn(scip, probdata, isFarkas, solfoundforday, vars[sortedvars[j].index], dualvalues,
                             newtours[addedCols], &newtourslength[addedCols], &newtoursobj[addedCols], allnames[addedCols]))
        {
            addedCols++;
        }
    }
    /* if we did not find any improvement for the tours that correspond to variables with value > 0
     * we check other variables */
    if(addedCols == 0 && SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == 1)
    {
        for(j = 0; j < nvars - nsortedvars; j++)
        {
            if(addedCols == MAX_HEURISTIC_TOURS || j > nvars/2)
            {
                break;
            }
            assert(vars[backupvars[j].index] != NULL);
            if( SCIPvarGetUbLocal(vars[backupvars[j].index]) < 0.5 )
            {
                continue;
            }
            if(investigateColumn(scip, probdata, isFarkas, solfoundforday, vars[backupvars[j].index], dualvalues,
                                 newtours[addedCols], &newtourslength[addedCols], &newtoursobj[addedCols], allnames[addedCols]))
            {
                addedCols++;
            }
        }
    }
    for(j = 0; j < addedCols; j++) // add each new column/tour to the model
    {
        solutionwindows = NULL;
        tmpobj = computeObjValue(scip, modeldata, &solutionwindows, &isfeasible, newtours[j], &duration, newtourslength[j], newtours[j][newtourslength[j]]);
        assert(tmpobj == newtoursobj[j]);
        assert(isfeasible);
        SCIP_CALL(SCIPcreateColumn(scip, probdata, allnames[j], FALSE, newtoursobj[j], newtours[j], newtourslength[j], duration, solutionwindows, newtours[j][newtourslength[j]]));
        SCIP_CALL( freeSolutionWindowArray(scip, solutionwindows, newtourslength[j]) );
    }
    for(j = 0; j < MAX_HEURISTIC_TOURS; j++)
    {
        SCIPfreeBlockMemoryArray(scip, &newtours[j], modeldata->nC);
    }

    SCIPfreeBlockMemoryArray(scip, &newtours, MAX_HEURISTIC_TOURS);
    SCIPfreeBlockMemoryArray(scip, &newtourslength, MAX_HEURISTIC_TOURS);
    SCIPfreeBlockMemoryArray(scip, &newtoursobj, MAX_HEURISTIC_TOURS);

    SCIPfreeBlockMemoryArray(scip, &solfoundforday, modeldata->nDays);
    SCIPfreeBlockMemoryArray(scip, &dualvalues, modeldata->nC - 1 + modeldata->nDays);
    SCIPfreeBlockMemoryArray(scip, &sortedvars, nsortedvars);
    SCIPfreeBlockMemoryArray(scip, &backupvars, nvars - nsortedvars);
    return SCIP_OKAY;
}


/** calls local search pricing heuristic */
SCIP_RETCODE heuristicPricing(
    SCIP*                scip,
    SCIP_Bool            isFarkas        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
    )
{
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);
    SCIP_CALL (localSearchPricing(scip, probdata, isFarkas));

    return SCIP_OKAY;
}