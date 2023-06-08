/**@file   initial_vrp.c
 * @brief  Methods for generating initial tours
 * @author Lukas Schürmann, University Bonn
 */

#include <string.h>

#include "probdata_vrp.h"
#include "tools_data.h"
#include "vardata_vrp.h"
#include "pricer_vrp.h"
#include "tools_vrp.h"
#include "initial_vrp.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/scip.h"

#include "time.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the vrp, all variables which are created, and all
 * constraints.
 */

/**@} */


/**@name Local methods
 *
 * @{
 */


/** @return Cumulative lengths of all time windows of node */
static
int getAvailableTimeOfCustomer(
        model_data*     modelData,
        int             node
){
    modelWindow* timewindow;
    int totaltime;

    assert(node >= 0 && node < modelData->nC-1);

    timewindow = modelData->timeWindows[node];
    assert(timewindow != NULL);

    totaltime = 0;
    while (timewindow != NULL)
    {
        totaltime += (timewindow->end_t - timewindow->start_t);
        timewindow = timewindow->next;
    }
    return totaltime;
}

/** improves and adds a set of tours to SCIP */
SCIP_RETCODE addToursToModel(
        SCIP*           scip,
        SCIP_PROBDATA*  probdata,
        int*            dayofnode,
        int**           tours,
        int*            tourlength,
        double*         tourobj,
        char*           algoName,
        SCIP_Bool       callback
){
    int i,j;
    int nDays, nC;
    int day;
    double obj;
    SCIP_Bool isfeasible;
    char name[SCIP_MAXSTRLEN];
    char strtmp[SCIP_MAXSTRLEN];
    int duration;
    int numUnused;
    double optionalCost = 0;

    model_data* modeldata;

    solutionWindow*** solutionwindows;

    modeldata = probdata->modeldata;
    nDays = modeldata->nDays;
    nC = modeldata->nC;

    SCIP_CALL(SCIPallocBlockMemoryArray(scip, &solutionwindows, nDays));

    /* count number of unvisited customers */
    numUnused = 0;
    for(j = 0; j < nC - 1; j++)
    {
        if(dayofnode[j] == -1)
        {
            numUnused++;
        }
    }
    if(numUnused)
    {
        /* tries to add unvisited customers */
        SCIP_CALL(addUnvisitedNodes(scip, modeldata, tours, tourlength, tourobj, dayofnode, TRUE));
        SCIP_CALL(addUnvisitedNodes(scip, modeldata, tours, tourlength, tourobj, dayofnode, FALSE));
        numUnused = 0;
        for(j = 0; j < nC - 1; j++)
        {
            if(dayofnode[j] == -1)
            {
                numUnused++;
            }
        }
        if(numUnused)
        {
            /* Create variables even if these do not induce a feasible solution */
            for (day = 0; day < nDays; day++)
            {
                if (tourlength[day] > 0)
                {
                    SCIP_CALL(rearrangeTour(scip, modeldata, tours[day], tourlength[day], &tourobj[day], day));
                    (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%2d: ", algoName, day);
                    for (i = 0; i < tourlength[day]; i++)
                    {
                        (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", tours[day][i]);
                        strcat(name, strtmp);
                    }
                    obj = computeObjValue(scip, modeldata, &(solutionwindows[day]), &isfeasible, tours[day], &duration, tourlength[day], day);
                    assert(duration <= WORKTIME_LIMIT);
                    assert(isfeasible);
                    SCIP_CALL( SCIPcreateColumn(scip, probdata, name, TRUE, obj, tours[day], tourlength[day], duration, solutionwindows[day], day));
                    SCIP_CALL( freeSolutionWindowArray(scip, solutionwindows[day], tourlength[day]) );
                }
            }
        }
    }
    if(!numUnused)
    {
        /* improve feasible solution */
        SCIP_CALL( twoNodeShift(scip, modeldata, tours, tourlength, tourobj, dayofnode));

        for (day = 0; day < nDays; day++)
        {
            /* create variables of feasible solution and add them to SCIP */
            if (tourlength[day] > 0)
            {
                (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "%s_%2d: ", algoName, day);
                for (i = 0; i < tourlength[day]; i++)
                {
                    (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", tours[day][i]);
                    strcat(name, strtmp);
                }
                obj = computeObjValue(scip, modeldata, &(solutionwindows[day]), &isfeasible, tours[day], &duration, tourlength[day], day);
                assert(duration <= WORKTIME_LIMIT);
                assert(obj == tourobj[day]);
                assert(isfeasible);
                SCIP_CALL( SCIPcreateColumn(scip, probdata, name, TRUE, tourobj[day], tours[day], tourlength[day], duration, solutionwindows[day], day));
                SCIP_CALL( freeSolutionWindowArray(scip, solutionwindows[day], tourlength[day]) );
            }
        }
    }
    SCIPfreeBlockMemoryArray(scip, &solutionwindows, nDays);

    /* check, if this incomplete solution is better than all others and update */
    for(j = 0; j < nC - 1; j++)
    {
        if(dayofnode[j] == -1)
        {
            optionalCost += modeldata->obj[j];
        }
    }
    if (probdata->useOptionals == FALSE && optionalCost < probdata->optionalCost)
    {
        probdata->optionalCost = optionalCost;
        for(j = 0; j < nC - 1; j++)
        {
            if(dayofnode[j] == -1)
            {
                probdata->optionalCustomers[j] = TRUE;
            } else {
                probdata->optionalCustomers[j] = FALSE;
            }
        }
    }
    return SCIP_OKAY;
}

/** Creates a tour by greedy decisions */
static
SCIP_RETCODE getGreedyTour(
        SCIP*           scip,
        model_data*     modeldata,
        int*            dayofnode,
        int*            tour,
        int             *length,
        int             day
){
    int tourlength;
    int currentNode;
    int gamma, i, nC;
    int bestTime;
    int arrivalTime;
    int latestTime;
    int* robusttimes;
    int* returntimes;
    int tourduration;
    int tourstart;
    SCIP_Bool isfeasible;
    modelWindow* timewindow;
    neighbor* nb;

    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &robusttimes, modeldata->maxDelayEvents + 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &returntimes, modeldata->maxDelayEvents + 1) );

    nC = modeldata->nC;
    /* start each tour at the depot */
    currentNode = nC - 1;
    tourlength = 0;
    tourstart = modeldata->shift_end;
    for (gamma = 0; gamma <= modeldata->maxDelayEvents; gamma++) /* initialize times at depot - no delays possible */
    {
        robusttimes[gamma] = modeldata->shift_start;
    }
    while (TRUE) /* extend tour while there are still some reachable customers */
    {
        bestTime = modeldata->shift_end;
        if(modeldata->neighbors != NULL){
            nb = modeldata->neighbors[currentNode][day];
            while (nb != NULL) /* iterate over the neighborhood of the head of the current tour */
            {
                i = nb->id;
                if (i == nC - 1 || dayofnode[i] >= 0) /* skip the depot and nodes that are already dayofnode by this or a tour of another day */
                {
                    nb = nb->next;
                    continue;
                }
                /* find the first time window of customer i that we can visit him in */
                arrivalTime = robusttimes[0] + modeldata->t_service[currentNode] +
                        getTravelTime(modeldata, currentNode, i, robusttimes[0] + modeldata->t_service[currentNode]);
                timewindow = getNextTimeWindow(modeldata, i, day, arrivalTime);
                if (timewindow == NULL) /* skip if there is none */
                {
                    nb = nb->next;
                    continue;
                }
                /* check at which time we would visit the customer if we got maximum delays (robustness)
                 * latestTime = -1 if it lays after the time window */
                latestTime = getLatestArrival(modeldata, timewindow, robusttimes, currentNode, i);

                if (latestTime == -1) /* skip in the latter case */
                {
                    nb = nb->next;
                    continue;
                }
                /* Greedy: best neighbor = neighbor with earliest worst-case arrvial time */
                for(gamma = 0; gamma <= modeldata->maxDelayEvents; gamma++)
                {
                    returntimes[gamma] = robusttimes[gamma];
                }
                updateRobustTimes(modeldata, returntimes, currentNode, i, day, &isfeasible, NULL, FALSE, NULL);
                timewindow = getNextTimeWindow(modeldata, nC - 1, day, returntimes[0]);
                assert(isfeasible);

                if (latestTime < bestTime && (getLatestArrival(modeldata, timewindow, returntimes, i, nC - 1) - tourstart < WORKTIME_LIMIT))
                {
                    bestTime = latestTime;
                    tour[tourlength] = i;
                }
                nb = nb->next;
            }
        }else{
            for(i = 0; i < nC - 1; i++)
            {
                if (i == currentNode || dayofnode[i] >= 0) /* skip the depot and nodes that are already dayofnode by this or a tour of another day */
                {
                    continue;
                }
                /* find the first time window of customer i that we can visit him in */
                arrivalTime = robusttimes[0] + modeldata->t_service[currentNode] +
                        getTravelTime(modeldata, currentNode, i, robusttimes[0] + modeldata->t_service[currentNode]);
                timewindow = getNextTimeWindow(modeldata, i, day, arrivalTime);
                if (timewindow == NULL) /* skip if there is none */
                {
                    continue;
                }
                /* check at which time we would visit the customer if we got maximum delays (robustness)
                 * latestTime = -1 if it lays after the time window */
                latestTime = getLatestArrival(modeldata, timewindow, robusttimes, currentNode, i);

                if (latestTime == -1) /* skip in the latter case */
                {
                    continue;
                }

                for(gamma = 0; gamma <= modeldata->maxDelayEvents; gamma++)
                {
                    returntimes[gamma] = robusttimes[gamma];
                }
                updateRobustTimes(modeldata, returntimes, currentNode, i, day, &isfeasible, NULL, FALSE, NULL);
                timewindow = getNextTimeWindow(modeldata, nC - 1, day, returntimes[0]);
                assert(isfeasible);
                /* Greedy: best neighbor = neighbor with earliest worst-case arrvial time */
                if (latestTime < bestTime && (getLatestArrival(modeldata, timewindow, returntimes, i, nC - 1) - tourstart < WORKTIME_LIMIT))
                {
                    bestTime = latestTime;
                    tour[tourlength] = i;
                }
            }
        }
        if (bestTime >= modeldata->shift_end) /* no valid neighbor found */
        {
            break;
        }
        /* extend the tour by best neighbor and update the robust arrival times */
        updateRobustTimes(modeldata, robusttimes, currentNode, tour[tourlength], day, &isfeasible, NULL, FALSE, NULL);
        assert(isfeasible);
        if(currentNode == nC - 1)
        {
            tourstart = robusttimes[0] - getTravelTime(modeldata, currentNode, tour[tourlength], robusttimes[0]);
        }
        currentNode = tour[tourlength];
        dayofnode[currentNode] = day;
        tourlength++;
    }
    *length = tourlength;

    if(tourlength > 0)
    {
        /* back to the depot */
        updateRobustTimes(modeldata, robusttimes, tour[tourlength - 1], nC - 1, day, &isfeasible, NULL, FALSE, NULL);
        tourduration = robusttimes[modeldata->maxDelayEvents] - tourstart;
        assert(tourduration >= 0);
        if(tourduration > WORKTIME_LIMIT || !isfeasible)
        {
            *length = 0;
            for(i = 0; i < tourlength; i++)
            {
                dayofnode[tour[i]] = -1;
            }
        }
    }
    SCIPfreeBlockMemoryArray(scip, &robusttimes, modeldata->maxDelayEvents + 1);
    SCIPfreeBlockMemoryArray(scip, &returntimes, modeldata->maxDelayEvents + 1);

    return SCIP_OKAY;
}


/** Adds a set of valid variables to the SPF - find tours in greedy manner */
SCIP_RETCODE createGreedyColumns(
        SCIP*                 scip,           /**< SCIP data structure */
        SCIP_PROBDATA*        probdata
)
{
   model_data* modeldata;
   int i, day, nC, nDays;

   int j,k;

   int** alltours;
   int* alltourlength;
   double* alltourobj;
   int* alltourduration;

   SCIP_Bool isfeasible;

   char algoName[] = "initGreedy";
   int* dayofnode;

    clock_t start, end;

   assert(probdata != NULL);
   assert(probdata->modeldata != NULL);

   modeldata = probdata->modeldata;

   nDays = modeldata->nDays;
   nC = modeldata->nC;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &dayofnode, nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alltourobj, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alltourlength, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alltours, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alltourduration, modeldata->nDays) );

    for(i = 0; i < nDays; i++)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &alltours[i], nC - 1) );
    }
    printf("Calculate and improve Greedy tours ... \n");

    start = clock();
    /* We want to find nDays vehicle routings. Therefore we start at a different day each time. */
    for (j = 0; j < nDays; j++)
    {
       for (i = 0; i < nC - 1; i++) /* We need to reset this status at each iteration */
       {
           dayofnode[i] = -1;
       }
       /* With the help of k, j and the modulo operation we get different permutations of the days
        * day(k) = (for j=0) 14, 13, ..., 0 -- (j=1)--> 0, 14, 13, ..., 1 -- (j=2)--> 1, 0, 14, 13, ..., 12 and so on */
       for (k = j + nDays - 1; k >= 0 + j; k--)
       {
           day = k % nDays;
           /* get tour in greedy manner */
           SCIP_CALL( getGreedyTour(scip, modeldata, dayofnode, alltours[day], &alltourlength[day], day));

           if(alltourlength[day] == 0)
           {
               alltourobj[day] = 0.0;
           }else{
               alltourobj[day] = computeObjValue(scip, modeldata, NULL, &isfeasible, alltours[day], NULL, alltourlength[day], day);
               assert(isfeasible);
           }
       }

        SCIP_CALL( addToursToModel( scip, probdata, dayofnode, alltours, alltourlength, alltourobj, algoName, FALSE));
        end = clock();
        if((double) (end - start) / CLOCKS_PER_SEC > MAX_GREEDY_TIME) // Do not invest too much time into the initial heuristics
        {
            printf("Max. greedy time exceeded. Continue ...\n");
            break;
        }
    }
    for(j=0;j < nDays; j++)
    {
        SCIPfreeBlockMemoryArray(scip, &alltours[j], nC - 1);
    }
    SCIPfreeBlockMemoryArray(scip, &alltours, nDays);
    SCIPfreeBlockMemoryArray(scip, &alltourlength, nDays);
    SCIPfreeBlockMemoryArray(scip, &alltourobj, nDays);
    SCIPfreeBlockMemoryArray(scip, &alltourduration, nDays);

   SCIPfreeBlockMemoryArray(scip, &dayofnode, nC - 1);
   return SCIP_OKAY;
}

/** dispatching heuristic to find feasible solution */
SCIP_RETCODE initialDispatching(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata
){
    int node, day, nC, nDays;
    int i, j, k;
    model_data* modelData;
    modelWindow* timewindow;
    inttuple* sortedCustomers;                  // customers sorted by total length of their time windows
    inttuple** sortedDaysForCustomer;           // days sorted by flexibility - for each customer
    int* sizeofDays;                            // value[day]: (number of available customers) - (number of customers already used for different days)
    int* nDaysOfCustomer;                       // value[customer]: number of potential days for customer
    int** tour;                                 // current tour for each day
    int* tourlength;                            // current tour length for each day
    double* tourobj;                            // current tour objective value for each day
    int unserved;                               // number of unserved customers;
    int* dayofnode;
    char algoName[] = "initDispatched";

    modelData = probdata->modeldata;
    nC = modelData->nC;
    nDays = modelData->nDays;

    /* Allocate memory */
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &dayofnode, nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sortedDaysForCustomer, nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &nDaysOfCustomer, nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sortedCustomers, nC - 1) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &sizeofDays, nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tour, nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tourlength, nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tourobj, nDays));

    /* Set up arrays */
    for(day = 0; day < nDays; day++)
    {
        tourlength[day] = 0;
        tourobj[day] = 0.0;
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(tour[day]), nC - 1));
        sizeofDays[day] = getNumberOfNeighbors(scip, modelData, nC - 1, day);
    }
    for(node = 0; node < nC - 1; node++)
    {
        dayofnode[node] = -1;
        sortedCustomers[node].index = node;
        sortedCustomers[node].value = getAvailableTimeOfCustomer(modelData, node);
        nDaysOfCustomer[node] = getNumActiveDaysOfCustomer(modelData, node);
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(sortedDaysForCustomer[node]), nDaysOfCustomer[node]) );
    }
    qsort(sortedCustomers, nC - 1, sizeof(sortedCustomers[0]), cmp_vrp);

    unserved = nC - 1;
    printf("Calculate and improve Dispatching tours ...\n");
    /* find a tour for each customer starting from the one with shortest time window length */
    for(i = nC - 2; i >= 0; i--)
    {
        node = sortedCustomers[i].index;    // current Customer

        /* Set up priority array for node */
        timewindow = modelData->timeWindows[node];
        for(day = 0; day < nDaysOfCustomer[node]; day++)
        {
            assert(timewindow != NULL);
            sortedDaysForCustomer[node][day].value = sizeofDays[timewindow->day];
            sortedDaysForCustomer[node][day].index = timewindow->day;
            timewindow = timewindow->next;
        }
        qsort(sortedDaysForCustomer[node], nDaysOfCustomer[node], sizeof(sortedDaysForCustomer[node][0]), cmp_vrp);

        /* Search for a day which tour can be extended by node starting at the day with the least flexibility */
        for(k = nDaysOfCustomer[node] - 1; k >= 0; k--)
        {
            day = sortedDaysForCustomer[node][k].index;     // current day

            /* check if tour[day] stays feasible when we add customer node and add if possible */
            if(addNodeToTour(scip, modelData, tour[day], &tourlength[day], &tourobj[day], node, day, NULL, SCIP_DEFAULT_INFINITY))
            {
                /* decrease flexibility of each day that node occurs in except for the current one */
                for(j = nDaysOfCustomer[node] - 1; j > k; j--)
                {
                    sizeofDays[sortedDaysForCustomer[node][j].index]--;
                }
                for(j = k-1; j >= 0; j--)
                {
                    sizeofDays[sortedDaysForCustomer[node][j].index]--;
                }
                unserved--;
                dayofnode[node] = day;
                break;
            }
            if(k == 0)
            {
                /* if node could not be added to any day decrease the flexibility of all potential days */
                for(j = 0; j < nDaysOfCustomer[node]; j++)
                {
                    sizeofDays[sortedDaysForCustomer[node][j].index]--;
                }
            }
        }
    }

    /* create and add variables for calculated tour */
    SCIP_CALL( addToursToModel( scip, probdata, dayofnode, tour, tourlength, tourobj, algoName, FALSE));

    /* free used memory */

    for(node = 0; node < nC - 1; node++)
    {
        SCIPfreeBlockMemoryArray(scip, &sortedDaysForCustomer[node], nDaysOfCustomer[node]);
    }
    for(day = 0; day < nDays; day++)
    {
        SCIPfreeBlockMemoryArray(scip, &tour[day], nC - 1);
    }
    SCIPfreeBlockMemoryArray(scip, &tourobj, nDays);
    SCIPfreeBlockMemoryArray(scip, &tourlength, nDays);
    SCIPfreeBlockMemoryArray(scip, &tour, nDays);
    SCIPfreeBlockMemoryArray(scip, &sizeofDays, nDays);
    SCIPfreeBlockMemoryArray(scip, &sortedCustomers, nC - 1);
    SCIPfreeBlockMemoryArray(scip, &nDaysOfCustomer, nC - 1);
    SCIPfreeBlockMemoryArray(scip, &sortedDaysForCustomer, nC - 1);
    SCIPfreeBlockMemoryArray(scip, &dayofnode, nC - 1);

    return SCIP_OKAY;
}

static
SCIP_RETCODE completeSolution(
        SCIP*                 scip,
        model_data*           modeldata,
        solution_data*        solutionData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode,
        SCIP_Bool*            complete
){
    int day, i;
    printf("Try to complete this partial solution... \n");
    
    *complete = TRUE;
    for(i = 0; i < modeldata->nC - 1; i++) // every new customer
    {
        if(dayofnode[i] == -1)
        {
            for(day = solutionData->nDays; day < modeldata->nDays; day++) // try to add this customer to one of the new days
            {
                if(addNodeToTour(scip, modeldata, tours[day], &tourlength[day], &tourobj[day], i, day, NULL, SCIP_DEFAULT_INFINITY))
                {
                    dayofnode[i] = day;
                    break;
                }
            }
            if(dayofnode[i] == -1)
            {
                for(day = 0; day < solutionData->nDays; day++) // try to add this customer to one of the old days
                {
                    if(addNodeToTour(scip, modeldata, tours[day], &tourlength[day], &tourobj[day], i, day, NULL, SCIP_DEFAULT_INFINITY))
                    {
                        dayofnode[i] = day;
                        break;
                    }
                }
                if(dayofnode[i] == -1)
                {
                    *complete = FALSE;
                }
            }
        }
    }
    
    if(!*complete) // if there is at least one customer that could not be added try greedy algorithm
    {
        *complete = TRUE;
        for(i = 0; i < modeldata->nC - 1; i++) // every new customer
        {
            if(dayofnode[i] >= solutionData->nDays)
            {
                dayofnode[i] = -1;
            }
        }
        for(day = solutionData->nDays; day < modeldata->nDays; day++) // call greedy algorithm for the new days
        {
            SCIP_CALL(getGreedyTour(scip, modeldata, dayofnode, tours[day], &tourlength[day], day));
        }
        for(i = 0; i < modeldata->nC - 1; i++) // every new customer
        {
            if(dayofnode[i] == -1)
            {
                *complete = FALSE;
            }
        }
    }
    if(!*complete) // if there is at least one customer that could not be added try to improve the partial solution
    {
        SCIP_CALL(addUnvisitedNodes(scip, modeldata, tours, tourlength, tourobj, dayofnode, TRUE));
        SCIP_CALL(addUnvisitedNodes(scip, modeldata, tours, tourlength, tourobj, dayofnode, FALSE));
    }
    *complete = TRUE;
    for(i = 0; i < modeldata->nC - 1; i++) // every new customer
    {
        if(dayofnode[i] == -1)
        {
            *complete = FALSE;
        }
    }
    return SCIP_OKAY;
}

/** creates and add variables for a given initial solution */
SCIP_RETCODE createSolutionColumns(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata,
        solution_data*        solutiondata
){
    int i, day;
    char algoName[] = "solutionInput";
    model_data* modeldata;
    int** tours;
    int* tourlength;
    double* tourobj;
    int* dayofnode;
    SCIP_Bool complete;
    SCIP_Bool isfeasible;
    int duration;
    double obj;
    modeldata = probdata->modeldata;

    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tours, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tourlength, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tourobj, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &dayofnode, modeldata->nC - 1) );

    for(i = 0; i < modeldata->nC - 1; i++)
    {
        dayofnode[i] = -1;
    }
    for(day = 0; day < solutiondata->nDays; day++) // extract the tours from solution data
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tours[day], modeldata->nC - 1) );
        tourlength[day] = solutiondata->tourLength[day];
        tourobj[day] = solutiondata->obj[day];

        for(i = 0; i < tourlength[day]; i++)
        {
            tours[day][i] = solutiondata->tour[day][i];
            dayofnode[solutiondata->tour[day][i]] = day;
        }
        obj = computeObjValue(scip, modeldata, NULL, &isfeasible, tours[day], &duration, tourlength[day], day);
        if (!isfeasible)
        {
            SCIPwarningMessage(scip, " *** Given initial solution for day %d is not feasible: ***\n", day);
            for (i = 0; i < tourlength[day]; i++)
            {
                printf(" %d - [%d, %d]\n", solutiondata->tour[day][i], solutiondata->timeWindows[solutiondata->tour[day][i]]->start_t, solutiondata->timeWindows[solutiondata->tour[day][i]]->latest_arrival);
            }
            printf("\n");
            printf("length = %d, given obj %f vs our obj %f\n", tourlength[day], tourobj[day], obj);
            assert(isfeasible);
        }
    }

    if(solutiondata->nDays == modeldata->nDays && solutiondata->nC == modeldata->nC) // if the solution is complete, just the corresponding columns
    {
        for(i = 0; i < modeldata->nC - 1; i++)
        {
            assert(dayofnode[i] >= 0);
        }
        SCIP_CALL( addToursToModel( scip, probdata, dayofnode, tours, tourlength, tourobj, algoName, FALSE));
    }
    else // else try to add missing/new customers
    {
        SCIP_Bool missingInBetween = FALSE;
        printf("Partial input solution detected!\n");
        assert(solutiondata->nC < modeldata->nC);
        for(i = 0; i < solutiondata->nC - 1; i++)
        {
            if (dayofnode[i] < 0)
            {
                missingInBetween = TRUE;
            }
        }
        for(i = solutiondata->nC - 1; i < modeldata->nC - 1; i++)
        {
            if (dayofnode[i] != 0)
            {
                missingInBetween = TRUE;
            }
        }
        if (missingInBetween)
        {
            printf("Detected missing customers in the middle of solution input.\n");
        }
        for(day = solutiondata->nDays; day < modeldata->nDays; day++)
        {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tours[day], modeldata->nC - 1) );
            tourlength[day] = 0;
            tourobj[day] = 0.0;
        }
        completeSolution(scip, modeldata, solutiondata, tours, tourlength, tourobj, dayofnode, &complete);
        if(complete)
        {
            printf("Success!\n");
        }
        else
        {
            printf("Failed completion! Trying heuristics without initial solution.\n");
            SCIP_CALL(createGreedyColumns(scip, probdata));
            SCIP_CALL(initialDispatching(scip, probdata));
        }
        SCIP_CALL( addToursToModel( scip, probdata, dayofnode, tours, tourlength, tourobj, algoName, FALSE));
    }

    return SCIP_OKAY;
}

/** Adds one empty route per day */
SCIP_RETCODE createInitialColumns(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_PROBDATA*        probdata            /**< problem data */
)
{
    char name[SCIP_MAXSTRLEN];
    int i, nC, nDays;
    int* tour;

    assert(probdata != NULL);
    assert(probdata->modeldata != NULL);

    nC = probdata->modeldata->nC;
    nDays = probdata->modeldata->nDays;
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tour, 1) );

    /* create variable with empty tour for each day */
    for( i = 0; i < nDays; i++ )
    {
        (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "initEmpty_%2d", i);
        tour[0] = nC - 1 + i;

        SCIP_CALL( SCIPcreateColumn(scip, probdata, name, TRUE, 0.0, tour, 0, 0, NULL, i));
    }
    SCIPfreeBlockMemoryArray(scip, &tour, 1);
    return SCIP_OKAY;
}
