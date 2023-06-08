/**@file   tools_vrp.c
 * @brief  debug parameters and diverse functions for column generation
 * @author Tim Niemann, TU Braunschweig
 */

#include <stdio.h>
#include <assert.h>

#include "tools_vrp.h"
#include "tools_data.h"
#include "probdata_vrp.h"
#include "postprocessing_vrp.h"
#include "pricer_vrp.h"


/**
 * compare function for qsort with a tuple structure
 */
int cmp_vrp(
    const void *a, 
    const void *b
    )
{
    tuple *a1 = (tuple *)a;
    tuple *a2 = (tuple *)b;
    if ((*a1).value > (*a2).value)
        return -1;
    else if ((*a1).value < (*a2).value)
        return 1;
    else
        return 0;
}

/**
 * Checks if the given node is the depot node 
 * @param modeldata underlying model data which contains the depot node
 * @param node node to check 
 * @return TRUE or FALSE */
SCIP_Bool nodeIsDepot(
    model_data*     modeldata,
    int             node
    )
{   
    assert(modeldata != NULL);
    assert(0 <= node && node < modeldata->nC);
    if (node == modeldata->nC - 1) {
        return TRUE;
    } else {
        return FALSE;
    }
}

/**
 * Computes the daytime dependet traveltime from one customer to the next
 * @param modeldata model data with all traveltimes
 * @param start start node
 * @param end end node
 * @param starttime starttime of travelling, in most cases endtime of service, this time is used to get the actual traveltime
 * @return daytime dependent traveltime from start to end */
int getTravelTime(
    model_data*     modeldata,
    int             start,
    int             end,
    int             starttime
    )
{
    int traveltime;

    assert(modeldata != NULL);
    assert(0 <= start && start < modeldata->nC);
    assert(0 <= end && end < modeldata->nC);
    return modeldata->t_travel[start][end];
    /* if time dependent travel times are deactivated, use the same value at every starttime */
    if (!TIME_DEPENDENT_TRAVEL_TIMES)
    {
        return modeldata->t_travel[start][end];
    }

    if (starttime < NOON_START)
    {
        assert(modeldata->t_travelAM != NULL);
        traveltime = modeldata->t_travelAM[start][end];
    } 
    else if (starttime < NOON_END)
    {
        assert(starttime >= NOON_START);
        assert(modeldata->t_travelNoon != NULL);
        traveltime = modeldata->t_travelNoon[start][end];
    } 
    else if (starttime < EVENING_START)
    {
        assert(starttime >= NOON_END);
        assert(modeldata->t_travelPM != NULL);
        traveltime = modeldata->t_travelPM[start][end];
    } else {
        assert(starttime >= EVENING_START);
        traveltime = modeldata->t_travel[start][end];
    }
    return traveltime;
}


/**
 * Find maximum between two or more integer variables
 * @param args Total number of integers
 * @param ... List of integer variables to find maximum
 * @return Maximum among all integers passed
 */
int max(
    int args, 
    ...
    )
{
    int i, max, cur;
    va_list valist;

    assert(args > 1);

    va_start(valist, args);
    max = INT_MIN;
    
    for(i=0; i<args; i++)
    {
        cur = va_arg(valist, int); // Get next elements in the list
        if(max < cur)
            max = cur;
    }
    
    va_end(valist); // Clean memory assigned by valist
    
    return max;
}

/** Gives the next possible time window for a service start at customer 'node' on day 'day' dependent on the current time
 * @param modeldata underlying model data
 * @param node given customer
 * @param day day for the time windows
 * @param arrivaltime arrival time without delays
 * @return time window of customer 'node' with lower bound = min { tw_lb | time window tw with tw_ub >= arrival time };
 *         NULL if node is the depot */
modelWindow* getNextTimeWindow(
        model_data*     modeldata,
        int             node,
        int             day,
        int             arrivaltime
    )
{
    modelWindow* timewindow;
    modelWindow* rttimewindow;

    assert(modeldata != NULL);
    assert(0 <= day && day < modeldata->nDays);
    assert(0 <= node && node < modeldata->nC);

    rttimewindow = NULL;

    timewindow = modeldata->timeWindows[node];
    while (timewindow != NULL)
    {
        if (timewindow->day != day)
        {
            timewindow = timewindow->next;
            continue;
        }
        assert(timewindow->day == day);
        /* lower bound of the actual time window */
        if (timewindow->start_t <= arrivaltime && timewindow->end_t >= arrivaltime)
        {
            rttimewindow = timewindow;
            break;
        }
        /* or lowest next time window, but search on for better ones */
        if(rttimewindow != NULL) {
            if (timewindow->start_t > arrivaltime && timewindow->start_t < rttimewindow->start_t)
            {
                rttimewindow = timewindow;
            }
        }else{
            if (timewindow->start_t > arrivaltime)
            {
                rttimewindow = timewindow;
            }
        }
        timewindow = timewindow->next;
    }
    return rttimewindow;
}

/** Gives the worst-case arrival time at node end,
 * when starting at node start and considering the current robust arrival times
 * @param modeldata underlying modeldata, contains the graph
 * @param timewindow the time window of customer end
 * @param robusttime robust arrival times of node start
 * @param start start node of the edge
 * @param end end node of the edge
 * @return if latest arrival time <= lower bound:       lower bound of time window
 *         else if latest arrival time <= upper bound:  latest arrival time
 *         else:                                        -1 */
int getLatestArrival(
        model_data*     modelData,
        modelWindow*    timewindow,
        int*            robusttime,
        int             start,
        int             end
    )
{
    int s,gamma,s_hat,t_hat;
    int no_dev, t_dev, s_dev, double_dev;
    int latestTime;
    int traveltime;

    assert(timewindow != NULL);

    s = modelData->t_service[start];
    s_hat = modelData->t_service_maxDev[start];
    t_hat = modelData->t_travel_maxDev[start][end];
    gamma = modelData->maxDelayEvents;
    traveltime = getTravelTime(modelData, start, end, robusttime[0] + s);

    if(gamma == 0)
    {
        latestTime = robusttime[gamma] + s + traveltime;
    }else if(gamma == 1)
    {
        no_dev = robusttime[gamma] + s + traveltime;
        s_dev = robusttime[gamma - 1] + s + traveltime + s_hat;
        t_dev = robusttime[gamma - 1] + s + traveltime + t_hat;
        latestTime = max(3, no_dev, t_dev, s_dev);
    }
    if(gamma >= 2)
    {
        no_dev = robusttime[gamma] + s + traveltime;
        s_dev = robusttime[gamma - 1] + s + traveltime + s_hat;
        t_dev = robusttime[gamma - 1] + s + traveltime + t_hat;
        double_dev = robusttime[gamma - 2] + s + traveltime + s_hat + t_hat;
        latestTime = max(4, no_dev, t_dev, s_dev, double_dev);
    }
    if(latestTime > timewindow->end_t)
    {
        return -1;
    }

    if(latestTime < timewindow->start_t)
    {
        return timewindow->start_t;
    }
    return latestTime;
}

/** @return Number of neighbors of node on a certain day */
int getNumberOfNeighbors(
        SCIP*           scip,
        model_data*     modelData,
        int             node,
        int             day
){
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;
    int numNeighbors;
    neighbor* nb;
    assert(node >= 0 && node < modelData->nC);
    assert(day >= 0 && day < modelData->nDays);

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);

    numNeighbors = 0;
    assert(modelData->neighbors != NULL);

    nb = pricerdata->neighbors[node][day];

    while (nb != NULL)
    {
        if(nb->id == modelData->nC - 1)
        {
            nb = nb->next;
            continue;
        }
        numNeighbors++;
        nb = nb->next;
    }
    return numNeighbors;
}

SCIP_RETCODE sortNeighborsOfNode(
        SCIP*           scip,
        model_data*     modelData,
        int             day,
        int             node,
        tuple*          list,
        int             nelements,
        double*         values
){
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;
    int i;
    neighbor* nb = NULL;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    /* Fill the sequence tupels */
    assert(modelData->neighbors != NULL);

    nb = pricerdata->neighbors[node][day];

    for (i = 0; i < nelements; i++)
    {
        if(nb == NULL)
        {
            printf("no neighbors for customer %d on day %d\n",node,day);
        }
        assert(nb != NULL);
        if (!nodeIsDepot(modelData, nb->id))
        {
            list[i].index = nb->id;
            list[i].value = values[nb->id];
        } else {
            i--;
        }
        nb = nb->next;
    }

    /* Sort the sequence by dualvalues */
    qsort(list, nelements, sizeof(list[0]), cmp_vrp);

    return SCIP_OKAY;
}

/** @return Number of days on which node could be served */
int getNumActiveDaysOfCustomer(
        model_data*     modelData,
        int             node
    )
{
    int* activeDays;
    modelWindow* timewindow;
    int numActiveDays;
    int day;

    assert(node >= 0 && node < modelData->nC);

    assert(modelData->timeWindows != NULL);

    activeDays = malloc(modelData->nDays * sizeof(int));
    numActiveDays = 0;

    for(day = 0; day < modelData->nDays; day++)
    {
        activeDays[day] = 0;
    }
    timewindow = modelData->timeWindows[node];
    while (timewindow != NULL)
    {
        activeDays[timewindow->day] = 1;
        timewindow = timewindow->next;
    }
    for(day = 0; day < modelData->nDays; day++)
    {
        if(activeDays[day])
        {
            numActiveDays++;
        }
    }
    free(activeDays);

    return numActiveDays;
}

/** Adds a node to a tour in the cheapest way.
 * @return returns FALSE if node could not be added */
SCIP_Bool addNodeToTour(
        SCIP*           scip,           /**< SCIP pointer */
        model_data*     modelData,      /**< underlying model data */
        int*            tour,           /**< tour (array of customers) that shall be investigated*/
        int             *tourlength,    /**< length of the tour */
        double          *tourobj,       /**< objective value of the tour */
        int             node,           /**< customer to be added */
        int             day,            /**< day of the tour */
        int             *newpos,        /**< to be set: position of the new customer in the tour */
        double          threshold       /**< only add node if additional costs lie under this value */
    )
{
    int i;
    int* tmptour;
    int length;
    SCIP_Bool isfeasible;
    double obj, bestobj;
    int bestpos;
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    assert(pricerdata != NULL);
    assert(pricerdata->isForbidden != NULL);

    assert(node >= 0 && node < modelData->nC);
    assert(day >= 0 && day < modelData->nDays);
    assert(*tourlength >= 0 && *tourlength < modelData->nC - 1);

    length = *tourlength;
    SCIP_CALL( SCIPallocMemoryArray(scip, &tmptour, modelData->nC ));
    bestpos = -1;
    bestobj = threshold > SCIP_DEFAULT_INFINITY / 2 ? threshold : *tourobj + threshold;

    /* initialize tour with new node at first position */
    tmptour[0] = node;
    for(i = 1; i < length + 1; i++)
    {
        tmptour[i] = tour[i - 1];
    }
    tmptour[length + 1] = modelData->nC - 1;
    if(!pricerdata->isForbidden[modelData->nC - 1][node] && !pricerdata->isForbidden[node][tmptour[1]])
    {
        obj = computeObjValue(scip, modelData, NULL, &isfeasible, tmptour, NULL, length + 1, day);
        if(isfeasible)
        {
            if(SCIPisSumPositive(scip, bestobj - obj))
            {
                bestpos = 0;
                bestobj = obj;
            }
        }
    }

    /* try every other position */
    for(i = 1; i < length + 1; i++){
        tmptour[i - 1] = tmptour[i];
        tmptour[i] = node;
        if(!pricerdata->isForbidden[tmptour[i - 1]][node] && !pricerdata->isForbidden[node][tmptour[i + 1]])
        {
            obj = computeObjValue(scip, modelData, NULL, &isfeasible, tmptour, NULL, length + 1, day);
            if(!isfeasible)
            {
                continue;
            }
            if(SCIPisSumPositive(scip, bestobj - obj))
            {
                bestpos = i;
                bestobj = obj;
            }
        }
    }
    SCIPfreeMemoryArray(scip, &tmptour);
    if (bestpos >= 0)
    {
        for(i = length; i > bestpos; i--)
        {
            tour[i] = tour[i-1];
        }
        tour[i] = node;
        (*tourlength)++;
        *tourobj = bestobj;
        if(newpos != NULL) *newpos = bestpos;
        return TRUE;
    }
    return FALSE;
}

SCIP_RETCODE rearrangeTour(
        SCIP*                 scip,
        model_data*           modelData,
        int*                  tour,
        int                   tourlength,
        double*               tourobj,
        int                   day
){
    int currentnode;
    int i, j;
    int pos;
    double obj;
    SCIP_Bool changed;
    SCIP_Bool isfeasible;
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    assert(pricerdata != NULL);
    assert(pricerdata->isForbidden != NULL);

    changed = TRUE;
    while(changed)
    {
        changed = FALSE;
        for(i = 0; i < tourlength; i++)
        {
            if(!isDeleteAllowed(tour, tourlength, modelData->nC - 1, i, pricerdata->isForbidden))
            {
                continue;
            }
            currentnode = tour[i];
            for(j = i; j > 0; j--)
            {
                tour[j] = tour[j - 1];
            }
            for(pos = 0; pos < tourlength; pos++)
            {
                if(isExchangeAllowed(tour, tourlength, modelData->nC - 1, pos, currentnode, pricerdata->isForbidden))
                {
                    tour[pos] = currentnode;
                    if(pos != i)
                    {
                        obj = computeObjValue(scip, modelData, NULL, &isfeasible, tour, NULL, tourlength, day);
                        if(isfeasible)
                        {
                            if(obj < *tourobj)
                            {
                                *tourobj = obj;
                                changed = TRUE;
                                break;
                            }
                        }
                    }
                }
                if(pos < tourlength - 1)
                {
                    tour[pos] = tour[pos + 1];
                }else{
                    for(j = pos; j > i; j--)
                    {
                        tour[j] = tour[j - 1];
                    }
                    tour[i] = currentnode;
                }
            }
        }
    }
    return SCIP_OKAY;
}

static
SCIP_RETCODE recursiveShiftNode(
        SCIP*               scip,
        model_data*         modelData,
        int**               tours,
        int*                tourlength,
        double*             tourobj,
        int*                dayofnode,
        int                 maxdepth,
        int                 currentdepth,
        double              sum,
        SCIP_Bool*          usedCustomers,
        SCIP_Bool*          usedDays,
        SCIP_Bool*          success
){
    int i, k;
    int customer;
    int oldday, newday;
    int oldpos, newpos;
    double oldobj_day1, newobj_day1;
    double oldobj_day2;
    double newsum;
    modelWindow* window;
    SCIP_Bool isfeasible;
    SCIP_Bool need_dependence, olddayUsed, newdayUsed;

    assert(currentdepth <= maxdepth);

    /* try out every customer to find an improvement */
    for(customer = 0; customer < modelData->nC - 1; customer++) //TODO: Alternatively, go through all tours -> no need to search for old pos
    {
        /* customer must not be used in prior recursions */
        if(!usedCustomers[customer])
        {
            usedCustomers[customer] = TRUE;
            oldday = dayofnode[customer];
            oldobj_day1 = tourobj[oldday];
            /* set the old day as used and remember if it already was before */
            olddayUsed = usedDays[oldday];
            usedDays[oldday] = TRUE;
            /* delete customer from current tour and remember its position */
            k = 0;
            tourlength[oldday]--;
            for(i = 0; i < tourlength[oldday]; i++)
            {
                if(tours[oldday][i] == customer)
                {
                    k = 1;
                    oldpos = i;
                }
                if(k == 1) tours[oldday][i] = tours[oldday][i + 1];
            }
            if(k == 0) oldpos = tourlength[oldday];
            /* calculate the new costs of this shorter tour */
            if(tourlength[oldday] == 0)
            {
                tourobj[oldday] = 0.0;
                isfeasible = TRUE;
            }
            else{
                tourobj[oldday] = computeObjValue(scip, modelData, NULL, &isfeasible, tours[oldday], NULL,
                                              tourlength[oldday], oldday);
            }
            newobj_day1 = tourobj[oldday]; // this needed if old day = new day
            if(!isfeasible)
            {
                /* undo everything */
                tourobj[oldday] = oldobj_day1;
                /* restore original tour of customer */
                for(i = tourlength[oldday]; i > oldpos; i--)
                {
                    tours[oldday][i] = tours[oldday][i - 1];
                }
                tours[oldday][oldpos] = customer;
                tourlength[oldday]++;
                /* set customer to unused */
                usedCustomers[customer] = FALSE;
                continue;
            }

            /* if current day of the customer is dependent of changed days, or it is the first recursion,
             * then the new day can be arbitrary
             * else the new day must depend on changed days */
            if(olddayUsed || currentdepth == 1)
            {
                need_dependence = FALSE;
            }else
            {
                need_dependence = TRUE;
            }
            /* try out every day the customer is available on */
            window = modelData->timeWindows[customer];
            while (window != NULL)
            {
                newday = window->day;
                /* if the old day was independent of the changed days the new must has to be dependent */
                if(need_dependence && !usedDays[newday])
                {
                    window = window->next;
                    continue;
                }
                /* if the final recursion has been reached we can make use of addNodeToTour to add customer
                 * to its cheapest position in the tour of new day.
                 * With the help of the threshold we can make sure that the customer only gets added,
                 * if the total objective change is negative */
                if(currentdepth == maxdepth)
                {
                    newsum = sum + newobj_day1 - oldobj_day1;
                    oldobj_day2 = tourobj[newday];
                    if(addNodeToTour(scip, modelData, tours[newday], &tourlength[newday], &tourobj[newday], customer,
                                     newday, NULL, -newsum))
                    {
                        *success = TRUE;
                        dayofnode[customer] = newday;
                        return SCIP_OKAY;
                    }
                /* else we check every possible position in this tour */
                }else{
                    /* set the new day as used and remember if it already was before */
                    newdayUsed = usedDays[newday];
                    usedDays[newday] = TRUE;

                    dayofnode[customer] = newday;
                    oldobj_day2 = tourobj[newday];
                    tourlength[newday]++;
                    /* try every spot of tour on new day - starting with position 0 */
                    for(i = tourlength[newday] - 1; i > 0; i--)
                    {
                        tours[newday][i] = tours[newday][i - 1];
                    }
                    tours[newday][0] = customer;
                    newpos = 0;
                    while (newpos <= tourlength[newday] - 1)
                    {
                        /* do not put customer back in the old spot */
                        if(newday != oldday || newpos != oldpos)
                        {
                            tourobj[newday] = computeObjValue(scip, modelData, NULL, &isfeasible, tours[newday], NULL,
                                                              tourlength[newday], newday);
                            if(isfeasible)
                            {
                                /* calculate new total change in costs */
                                newsum = sum + (newobj_day1 - oldobj_day1) + (tourobj[newday] - oldobj_day2);

                                if(currentdepth < maxdepth)
                                {
                                    assert(currentdepth < maxdepth);
                                    /* Dive deeper */
                                    SCIP_CALL(recursiveShiftNode(scip, modelData, tours, tourlength, tourobj, dayofnode,
                                                                 maxdepth, currentdepth + 1, newsum, usedCustomers, usedDays, success));
                                    if(*success) /* Improvement has been found */
                                    {
                                        return SCIP_OKAY;
                                    }
                                }else
                                {
                                    if(SCIPisSumNegative(scip, newsum))
                                    {
                                        *success = TRUE;
                                        return SCIP_OKAY;
                                    }
                                }

                            }
                        }
                        /* no improvement on this new position -> check next */
                        tours[newday][newpos] = tours[newday][newpos + 1];
                        newpos++;
                        tours[newday][newpos] = customer;
                    }
                    /* no improvement on this day -> check next */
                    tourobj[newday] = oldobj_day2;
                    tourlength[newday]--;
                    if(!newdayUsed) usedDays[newday] = FALSE;
                }
                window = window->next;
            }

            /* no improvement with this customer -> check next */
            tourobj[oldday] = oldobj_day1;
            /* restore original tour of customer */
            for(i = tourlength[oldday]; i > oldpos; i--)
            {
                tours[oldday][i] = tours[oldday][i - 1];
            }
            tours[oldday][oldpos] = customer;
            tourlength[oldday]++;
            /* set customer to unused */
            dayofnode[customer] = oldday;
            usedCustomers[customer] = FALSE;
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE recursiveShiftNode2(
        SCIP*               scip,
        model_data*         modelData,
        int**               tours,
        int*                tourlength,
        double*             tourobj,
        int*                dayofnode,
        int                 maxdepth,
        int                 currentdepth,
        double              sum,
        SCIP_Bool*          usedCustomers,
        SCIP_Bool*          usedDays,
        SCIP_Bool*          success
){
    int i, j, k;
    int customer;
    int oldday, newday;
    int oldpos, newpos;
    double oldobj_day1, newobj_day1;
    double oldobj_day2;
    double newsum;
    modelWindow* window;
    SCIP_Bool isfeasible;
    SCIP_Bool need_dependence, olddayUsed, newdayUsed;

    assert(currentdepth <= maxdepth);

    /* check every customer by iteration over all tours */
    for(j = 0; j < modelData->nDays; j++)
    {
        for(k = 0; k < tourlength[j]; k++)
        {
            customer = tours[j][k];
            if(!usedCustomers[customer])
            {
                usedCustomers[customer] = TRUE;
                /* save customer data */
                oldday = dayofnode[customer];
                oldobj_day1 = tourobj[oldday];
                oldpos = k;
                /* set the old day as used and remember if it already was before */
                olddayUsed = usedDays[oldday];
                usedDays[oldday] = TRUE;
                /* delete customer from current tour and remember its position */
                tourlength[oldday]--;
                for(i = k; i < tourlength[oldday]; i++)
                {
                    tours[oldday][i] = tours[oldday][i + 1];
                }
                /* calculate the new costs of this shorter tour */
                if(tourlength[oldday] == 0)
                {
                    tourobj[oldday] = 0.0;
                    isfeasible = TRUE;
                }
                else{
                    tourobj[oldday] = computeObjValue(scip, modelData, NULL, &isfeasible, tours[oldday], NULL,
                                                      tourlength[oldday], oldday);
                }
                newobj_day1 = tourobj[oldday]; // this needed if old day = new day
                if(!isfeasible)
                {
                    /* undo everything */
                    tourobj[oldday] = oldobj_day1;
                    /* restore original tour of customer */
                    for(i = tourlength[oldday]; i > oldpos; i--)
                    {
                        tours[oldday][i] = tours[oldday][i - 1];
                    }
                    tours[oldday][oldpos] = customer;
                    tourlength[oldday]++;
                    /* set customer to unused */
                    usedCustomers[customer] = FALSE;
                    continue;
                }

                /* if current day of the customer is dependent of changed days, or it is the first recursion,
                 * then the new day can be arbitrary
                 * else the new day must depend on changed days */
                if(olddayUsed || currentdepth == 1)
                {
                    need_dependence = FALSE;
                }else
                {
                    need_dependence = TRUE;
                }
                /* try out every day the customer is available on */
                window = modelData->timeWindows[customer];
                while (window != NULL)
                {
                    newday = window->day;
                    /* if the old day was independent of the changed days the new must has to be dependent */
                    if(need_dependence && !usedDays[newday])
                    {
                        window = window->next;
                        continue;
                    }
                    /* if the final recursion has been reached we can make use of addNodeToTour to add customer
                     * to its cheapest position in the tour of new day.
                     * With the help of the threshold we can make sure that the customer only gets added,
                     * if the total objective change is negative */
                    if(currentdepth == maxdepth)
                    {
                        newsum = sum + newobj_day1 - oldobj_day1;
                        oldobj_day2 = tourobj[newday];
                        if(addNodeToTour(scip, modelData, tours[newday], &tourlength[newday], &tourobj[newday], customer,
                                         newday, NULL, -newsum))
                        {
                            *success = TRUE;
                            dayofnode[customer] = newday;
                            return SCIP_OKAY;
                        }
                        /* else we check every possible position in this tour */
                    }else{
                        /* set the new day as used and remember if it already was before */
                        newdayUsed = usedDays[newday];
                        usedDays[newday] = TRUE;

                        dayofnode[customer] = newday;
                        oldobj_day2 = tourobj[newday];
                        tourlength[newday]++;
                        /* try every spot of tour on new day - starting with position 0 */
                        for(i = tourlength[newday] - 1; i > 0; i--)
                        {
                            tours[newday][i] = tours[newday][i - 1];
                        }
                        tours[newday][0] = customer;
                        newpos = 0;
                        while (newpos <= tourlength[newday] - 1)
                        {
                            /* do not put customer back in the old spot */
                            if(newday != oldday || newpos != oldpos)
                            {
                                tourobj[newday] = computeObjValue(scip, modelData, NULL, &isfeasible, tours[newday], NULL,
                                                                  tourlength[newday], newday);
                                if(isfeasible)
                                {
                                    /* calculate new total change in costs */
                                    newsum = sum + (newobj_day1 - oldobj_day1) + (tourobj[newday] - oldobj_day2);

                                    assert(currentdepth < maxdepth);
                                    /* Dive deeper */
                                    SCIP_CALL(recursiveShiftNode2(scip, modelData, tours, tourlength, tourobj, dayofnode,
                                                                 maxdepth, currentdepth + 1, newsum, usedCustomers, usedDays, success));
                                    if(*success) /* Improvement has been found */
                                    {
                                        return SCIP_OKAY;
                                    }
                                }
                            }
                            /* no improvement on this new position -> check next */
                            tours[newday][newpos] = tours[newday][newpos + 1];
                            newpos++;
                            tours[newday][newpos] = customer;
                        }
                        /* no improvement on this day -> check next */
                        tourobj[newday] = oldobj_day2;
                        tourlength[newday]--;
                        if(!newdayUsed) usedDays[newday] = FALSE;
                    }
                    window = window->next;
                }

                /* no improvement with this customer -> check next */
                tourobj[oldday] = oldobj_day1;
                /* restore original tour of customer */
                for(i = tourlength[oldday]; i > oldpos; i--)
                {
                    tours[oldday][i] = tours[oldday][i - 1];
                }
                tours[oldday][oldpos] = customer;
                tourlength[oldday]++;
                /* set customer to unused */
                dayofnode[customer] = oldday;
                usedCustomers[customer] = FALSE;
            }
        }
    }

    return SCIP_OKAY;
}

SCIP_RETCODE localSearch(
        SCIP*               scip,
        model_data*         modelData,
        int**               tours,
        int*                tourlength,
        double*             tourobj,
        int*                dayofnode,
        int                 maxdepth
){
    int i, k;
    SCIP_Bool* usedCustomers;           /* We do not want move one customer multiple times in a single search */
    SCIP_Bool* usedDays;                /* Each movement needs to correlate with at least one prior used day */
    SCIP_Bool success;

    SCIP_CALL( SCIPallocMemoryArray(scip, &usedCustomers, modelData->nC - 1) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &usedDays, modelData->nDays) );

    k = 1;
    while(k <= maxdepth)
    {
        /* initialization of a new search */
        for(i = 0; i < modelData->nC - 1; i++)
        {
            usedCustomers[i] = FALSE;
        }
        for(i = 0; i < modelData->nDays; i++)
        {
            usedDays[i] = FALSE;
        }
        success = FALSE;

        SCIP_CALL(recursiveShiftNode(scip, modelData, tours, tourlength, tourobj, dayofnode,
                                     k, 1, 0.0, usedCustomers, usedDays, &success));
        /* reset the recursive depth to 1 if an improvement has been found
         * else increase it to expand the search */
        if(success)
        {
            k = 1;
        }else
        {
            k++;
        }
    }

    SCIPfreeMemoryArray(scip, &usedCustomers);
    SCIPfreeMemoryArray(scip, &usedDays);

    return SCIP_OKAY;
}

SCIP_Bool isDeleteAllowed(
        int*                tour,
        int                 length,
        int                 depot,
        int                 pos,
        SCIP_Bool**         isForbidden
)
{
    assert(length >= 1);
    if(pos == 0)
    {
        if(length == 1) return TRUE;
        if(isForbidden[depot][tour[pos + 1]]) return FALSE;
    }
    else if(pos == length - 1)
    {
        if(isForbidden[tour[pos - 1]][depot]) return FALSE;
    }
    else
    {
        if(isForbidden[tour[pos - 1]][tour[pos + 1]]) return FALSE;
    }
    return TRUE;
}

SCIP_Bool isExchangeAllowed(
        int*                tour,
        int                 length,
        int                 depot,
        int                 pos,
        int                 new,
        SCIP_Bool**         isForbidden
)
{
    assert(length >= 1);
    if(pos == 0)
    {
        if(isForbidden[depot][new]) return FALSE;
        if(length == 1 && isForbidden[new][depot]) return FALSE;
        else if(length > 1 && isForbidden[new][tour[1]]) return FALSE;
    }
    else if(pos == length - 1)
    {
        if(isForbidden[new][depot]) return FALSE;
        if(isForbidden[tour[pos - 1]][new]) return FALSE;
    }
    else
    {
        if(isForbidden[tour[pos - 1]][new]) return FALSE;
        if(isForbidden[new][tour[pos + 1]]) return FALSE;
    }
    return TRUE;
}

/**
 * As long as there are customers u, v, w such that exchanging v by u, exchanging w by v
 * and giving w a new position in one of the tours yields a better solution,
 * apply this step.
 * */
SCIP_RETCODE twoNodeShift(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode
)
{
    double bestimprovement;
    int currentnode, node, pos;
    int nC;
    int tmplength;
    int* tourwithoutnode;
    int* tmptour;
    int k;
    SCIP_Bool isfeasible;
    SCIP_Bool changed;
    int day1, day2, day3;
    int pos1, newpos;
    double oldobj1, oldobj2, oldobj3;
    double newobj1, newobj2, newobj3;
    double oldtotal, newtotal;

    double bestnewobj1, bestnewobj2;
    int bestnewday1, bestnewday2;
    int bestnewpos1, bestnewpos2;
    int numimprov = 0;

    modelWindow* window;
    modelWindow* window2;
    SCIP_PRICER* pricer;
    SCIP_PRICERDATA* pricerdata;

    pricer = SCIPfindPricer(scip, "vrp");
    assert(pricer != NULL);
    pricerdata = SCIPpricerGetData(pricer);
    assert(pricerdata != NULL);
    assert(pricerdata->isForbidden != NULL);

    assert(modelData != NULL);
    assert(modelData->timeWindows != NULL);
    nC = modelData->nC;

    SCIP_CALL(SCIPallocBlockMemoryArray(scip, &tourwithoutnode, nC - 1) );
    SCIP_CALL(SCIPallocBlockMemoryArray(scip, &tmptour, nC - 1) );

    changed = TRUE;
    while(changed){
        changed = FALSE;
        for (currentnode = 0; currentnode < nC - 1; currentnode++)
        {
            bestimprovement = 0.0;
            day1 = dayofnode[currentnode];
            assert(day1 >= 0);
            k = 0;
            oldobj1 = tourobj[day1];
            for (pos = 0; pos < tourlength[day1]; pos++)
            {
                if (tours[day1][pos] == currentnode)
                {
                    if(!isDeleteAllowed(tours[day1], tourlength[day1], modelData->nC - 1, pos, pricerdata->isForbidden))
                    {
                        k = -1;
                        break;
                    }
                    k = 1;
                    continue;
                }
                tourwithoutnode[pos - k] = tours[day1][pos];
            }
            if(k == -1)
            {
                continue;
            }
            assert(k == 1);
            tourlength[day1]--;
            if (tourlength[day1] == 0)
            {
                newobj1 = 0.0;
                isfeasible = TRUE;
            } else {
                newobj1 = computeObjValue(scip, modelData, NULL, &isfeasible, tourwithoutnode, NULL, tourlength[day1], day1);
            }
            if(!isfeasible)
            {
                tourlength[day1]++;
                continue;
            }

            window = modelData->timeWindows[currentnode];
            assert(window != NULL);
            while (window != NULL)
            {
                day2 = window->day;
                for (pos1 = 0; pos1 < tourlength[day2]; pos1++) {
                    oldobj2 = tourobj[day2];
                    if(day2 == day1)
                    {
                        if(!isExchangeAllowed(tourwithoutnode, tourlength[day2], modelData->nC - 1, pos1, currentnode,
                                              pricerdata->isForbidden))
                        {
                            continue;
                        }
                        node = tourwithoutnode[pos1];
                        tourwithoutnode[pos1] = currentnode;

                        newobj2 = computeObjValue(scip, modelData, NULL, &isfeasible, tourwithoutnode, NULL, tourlength[day2],
                                                  day2);
                        if (!isfeasible)
                        {
                            tourwithoutnode[pos1] = node;
                            continue;
                        }
                    }else {
                        if(!isExchangeAllowed(tours[day2], tourlength[day2], modelData->nC - 1, pos1, currentnode,
                                              pricerdata->isForbidden))
                        {
                            continue;
                        }
                        node = tours[day2][pos1];
                        tours[day2][pos1] = currentnode;
                        newobj2 = computeObjValue(scip, modelData, NULL, &isfeasible, tours[day2], NULL, tourlength[day2], day2);
                        if (!isfeasible)
                        {
                            tours[day2][pos1] = node;
                            continue;
                        }
                    }
                    window2 = modelData->timeWindows[node];
                    while (window2 != NULL) {
                        day3 = window2->day;
                        for(pos = 0; pos < tourlength[day3]; pos++)
                        {
                            if(day3 == day1)
                            {
                                tmptour[pos] = tourwithoutnode[pos];
                            }else{
                                tmptour[pos] = tours[day3][pos];
                            }
                        }
                        tmplength = tourlength[day3];
                        oldobj3 = tourobj[day3];
                        if (addNodeToTour(scip, modelData, tmptour, &tmplength, &newobj3, node, day3, &newpos, SCIP_DEFAULT_INFINITY))
                        {
                            if(day1 == day2)
                            {
                                if(day2 == day3)
                                {
                                    oldtotal = oldobj1;
                                    newtotal = newobj3;
                                }else{
                                    oldtotal = oldobj1 + oldobj3;
                                    newtotal = newobj2 + newobj3;
                                }
                            }else if(day1 == day3)
                            {
                                oldtotal = oldobj1 + oldobj2;
                                newtotal = newobj2 + newobj3;
                            }else if(day2 == day3)
                            {
                                oldtotal = oldobj1 + oldobj2;
                                newtotal = newobj1 + newobj3;
                            }else{
                                oldtotal = oldobj1 + oldobj2 + oldobj3;
                                newtotal = newobj1 + newobj2 + newobj3;
                            }
                            if(oldtotal - newtotal > bestimprovement)
                            {
                                bestimprovement = oldtotal - newtotal;
                                bestnewpos1 = pos1;
                                bestnewpos2 = newpos;
                                bestnewday1 = day2;
                                bestnewday2 = day3;
                                bestnewobj1 = newobj2;
                                bestnewobj2 = newobj3;
                            }
                        }
                        window2 = window2->next;
                    }
                    if(day1 == day2)
                    {
                        tourwithoutnode[pos1] = node;
                    }else{
                        tours[day2][pos1] = node;
                    }
                }
                window = window->next;
            }
            if(bestimprovement > 0.0)
            {
                changed = TRUE;
                for(pos = 0; pos < tourlength[day1]; pos++)
                {
                    tours[day1][pos] = tourwithoutnode[pos];
                }
                tourobj[day1] = newobj1;
                node = tours[bestnewday1][bestnewpos1];
                tours[bestnewday1][bestnewpos1] = currentnode;
                tourobj[bestnewday1] = bestnewobj1;
                for(pos = tourlength[bestnewday2]; pos > bestnewpos2; pos--)
                {
                    tours[bestnewday2][pos] = tours[bestnewday2][pos - 1];
                }
                tours[bestnewday2][bestnewpos2] = node;
                tourobj[bestnewday2] = bestnewobj2;
                tourlength[bestnewday2]++;
                dayofnode[currentnode] = bestnewday1;
                dayofnode[node] = bestnewday2;
                numimprov++;

            }else{
                tourlength[day1]++;
            }
        }
    }
    SCIPfreeBlockMemoryArray(scip, &tourwithoutnode, nC - 1);
    SCIPfreeBlockMemoryArray(scip, &tmptour, nC - 1);
    return SCIP_OKAY;
}

/**
 * As long as there are customers u, v such that exchanging v by u and giving v a new position in one of the tours
 * yields a better solution, apply this step.
 * */
SCIP_RETCODE shiftNodes(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode
    )
{
    double oldobj;
    double newobj;
    double currentobj;
    double newtourobj;
    double objchange;
    double bestimprovement;
    int currentnode, node, day, pos;
    int nC, nDays;
    int tmplength;
    int bestDay;
    int* tmptour;
    int k, newpos, bestpos;
    double bestobj;
    SCIP_Bool isfeasible;
    SCIP_Bool changed;
    int* visited;

    modelWindow* window;

    assert(modelData != NULL);
    assert(modelData->timeWindows != NULL);

    nC = modelData->nC;
    nDays = modelData->nDays;

    SCIP_CALL(SCIPallocMemoryArray(scip, &tmptour, nC - 1) );
    SCIP_CALL(SCIPallocMemoryArray(scip, &visited, nC - 1) );

    currentobj = 0.0;
    for(day = 0; day < nDays; day++)
    {
        currentobj += tourobj[day];
    }
    changed = TRUE;
    while(changed)
    {
        changed = FALSE;
        for (currentnode = 0; currentnode < nC - 1; currentnode++)
        {
            visited[currentnode] = 0;
            bestimprovement = 0.0;
            day = dayofnode[currentnode];
            assert(day >= 0);
            k = 0;
            tmplength = tourlength[day] - 1;
            newtourobj = tourobj[day];

            for (pos = 0; pos < tourlength[day]; pos++)
            {
                if (tours[day][pos] == currentnode)
                {
                    k = 1;
                    continue;
                }
                tmptour[pos - k] = tours[day][pos];
            }
            if (tmplength == 0)
            {
                objchange = tourobj[day];
                newtourobj = 0.0;
                isfeasible = TRUE;
            } else {
                newtourobj = computeObjValue(scip, modelData, NULL, &isfeasible, tmptour, NULL, tmplength, day);
                objchange = tourobj[day] - newtourobj;
            }
            if(!isfeasible)
            {
                continue;
            }

            window = modelData->timeWindows[currentnode];
            assert(window != NULL);
            while (window != NULL)
            {
                day = window->day;
                oldobj = tourobj[day];
                tmplength = tourlength[day];
                k = 0;
                for (pos = 0; pos < tourlength[day]; pos++)
                {
                    if (tours[day][pos] == currentnode)
                    {
                        k = 1;
                        oldobj = newtourobj;
                        tmplength--;
                        continue;
                    }
                    tmptour[pos - k] = tours[day][pos];
                }
                if (addNodeToTour(scip, modelData, tmptour, &tmplength, &newobj, currentnode, day, &newpos, SCIP_DEFAULT_INFINITY))
                {
                    if (objchange - (newobj - oldobj) > bestimprovement)
                    {
                        bestimprovement = objchange - (newobj - oldobj);
                        bestobj = newobj;
                        bestDay = day;
                        bestpos = newpos;
                    }
                }
                window = window->next;
            }
            if (bestimprovement > 0.0)
            {
                changed = TRUE;
                day = dayofnode[currentnode];
                k = 0;
                for (pos = 0; pos < tourlength[day] - 1; pos++)
                {
                    if (k == 1)
                    {
                        tours[day][pos] = tours[day][pos + 1];
                    } else if (tours[day][pos] == currentnode)
                    {
                        k = 1;
                        tours[day][pos] = tours[day][pos + 1];
                        continue;
                    }
                }
                tourobj[day] = newtourobj;
                tourlength[day]--;
                dayofnode[currentnode] = bestDay;

                for (pos = tourlength[bestDay]; pos > bestpos; pos--)
                {
                    tours[bestDay][pos] = tours[bestDay][pos - 1];
                }
                tours[bestDay][bestpos] = currentnode;
                tourlength[bestDay]++;
                tourobj[bestDay] = bestobj;
            }
        }
    }
//    printf("Obj Before: %f\n", currentobj);
    currentobj = 0.0;
    for(day = 0; day < nDays; day++)
    {
        currentobj += tourobj[day];
        for(pos = 0; pos < tourlength[day]; pos++)
        {
            visited[tours[day][pos]]++;
        }
    }
    for(node = 0; node < nC - 1; node++)
    {
        if(visited[node] == 0)
        {
            assert(FALSE);
        }else if(visited[node] > 1)
        {
            assert(FALSE);
        }
    }

    SCIPfreeMemoryArray(scip, &tmptour);
    return SCIP_OKAY;
}

/** Tries to add unvisited customers to a incomplete vehicle routing by moving assigned customers to different spots */
extern
SCIP_RETCODE addUnvisitedNodes(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode,
        SCIP_Bool             isGreedy
    )
{
    int node, currentnode, day, nC;
    int day2;
    int pos;
    double obj;
    modelWindow* window;
    modelWindow* window2;
    SCIP_Bool isfeasible;
    SCIP_Bool foundExchange;
    int duration;

    assert(modelData != NULL);
    assert(modelData->timeWindows != NULL);
    nC = modelData->nC;

    for(currentnode = 0; currentnode < nC - 1; currentnode++)
    {
        if(dayofnode[currentnode] == -1) { // check every unvisited customer
            window = modelData->timeWindows[currentnode];
            foundExchange = FALSE;
            while (window != NULL) {
                day = window->day;
                if(isGreedy) // first: try out every possible spot in the daily tour - already checked in dispatching heuristic
                {
                    if (addNodeToTour(scip, modelData, tours[day], &tourlength[day], &tourobj[day], currentnode, day,
                                      NULL, SCIP_DEFAULT_INFINITY))
                    {
                        dayofnode[currentnode] = day;
                        break;
                    }
                }
                for (pos = 0; pos < tourlength[day]; pos++) // try to exchange customer
                {
                    node = tours[day][pos];
                    tours[day][pos] = currentnode;
                    obj = computeObjValue(scip, modelData, NULL, &isfeasible, tours[day], &duration, tourlength[day], day);
                    if (!isfeasible)
                    {
                        tours[day][pos] = node;
                        continue;
                    }
                    window2 = modelData->timeWindows[node];
                    while (window2 != NULL) // find a new spot for the exchanged customer
                    {
                        day2 = window2->day;
                        if (addNodeToTour(scip, modelData, tours[day2], &tourlength[day2], &tourobj[day2], node,
                                          day2, NULL, SCIP_DEFAULT_INFINITY))
                        { // if a valid new spot has been found, apply this exchange
                            foundExchange = TRUE;
                            dayofnode[currentnode] = day;
                            dayofnode[node] = day2;
                            if(day != day2)
                            {
                                tourobj[day] = obj;
                            }
                            break;
                        }
                        window2 = window2->next;
                    }
                    if (foundExchange) break;
                    tours[day][pos] = node; // undo exchange if no new valid spot has been found
                }

                if(foundExchange) break;
                window = window->next;
            }
        }

    }


    return SCIP_OKAY;
}

/** Converts the robust arrival times of customer start
 * into the robust arrival times of customer end by using the arc (start, end) */
SCIP_RETCODE updateRobustTimes(
    model_data* modelData,
    int*        robustTimes,
    int         start,
    int         end,
    int         day,
    SCIP_Bool   *isfeasible,
    double      *windowweight,
    SCIP_Bool   gethardtw,
    modelWindow** hardwindow
    )
{
    int gamma, earliestService, serviceAndTraveltime;
    int serviceEnd = robustTimes[0] + modelData->t_service[start];
    int traveltime = getTravelTime(modelData, start, end, serviceEnd);
    int earliestArrival = serviceEnd + traveltime;
    int no_dev, t_dev, s_dev, double_dev;
    modelWindow* window = getNextTimeWindow(modelData, end, day, earliestArrival);

    *isfeasible = TRUE;
    assert(modelData != NULL);
    assert(robustTimes != NULL);

    if (window == NULL)
    {
        *isfeasible = FALSE;
        return SCIP_OKAY;
    }
    if (windowweight != NULL)
    {
        *windowweight = window->weigth;
    }
    if(gethardtw)
    {
        *hardwindow = window;
    }
    earliestService = window->start_t;
    serviceAndTraveltime = modelData->t_service[start] + traveltime;
    no_dev = serviceAndTraveltime;
    t_dev = serviceAndTraveltime + modelData->t_travel_maxDev[start][end];
    s_dev = serviceAndTraveltime + modelData->t_service_maxDev[start];
    double_dev = serviceAndTraveltime + modelData->t_service_maxDev[start] + modelData->t_travel_maxDev[start][end];
    for (gamma = modelData->maxDelayEvents; gamma >= 0; gamma--)
    {
        /* Case: startnode is the depot */
        if (start == modelData->nC - 1)
        {
            robustTimes[gamma] = max(2, earliestService, modelData->shift_start + traveltime);
            /* directly after the depot could only be a exactly one delay: travel time */
            if (gamma > 0)
            {
                robustTimes[gamma] += modelData->t_travel_maxDev[start][end];
            }
        }
        /* Case: arrivaltimes without deviations */
        else if (gamma == 0)
        {
            /* a_j = Y_{i0} + s_i + t_{ij} */
            robustTimes[gamma] = max(2, earliestService, earliestArrival);
        }
        /* Case: arrivaltimes with exactly one deviations */
        else if (gamma == 1)
        {
            robustTimes[gamma] = max(4, earliestService, no_dev + robustTimes[1], t_dev + robustTimes[0], s_dev + robustTimes[0]);
        }
        /* All other cases */
        else
        {
            robustTimes[gamma] = max(5, earliestService, no_dev + robustTimes[gamma], t_dev + robustTimes[gamma - 1],
                                     s_dev + robustTimes[gamma - 1], double_dev + robustTimes [gamma - 2]);
        }
    }

    /* check that arrival times are in ascending order */
    #ifndef NDEBUG
    for (gamma = 1; gamma <= modelData->maxDelayEvents; gamma++)
    {
        assert(robustTimes[gamma - 1] <= robustTimes[gamma]);
    }
    #endif

    /* check that arrival times are all inside a feasible time window */
    if (window->start_t > robustTimes[0] || window->end_t < robustTimes[modelData->maxDelayEvents])
    {
        *isfeasible = FALSE;
    }

    return SCIP_OKAY;
}

double computeObjValue(
        SCIP*                 scip,
        model_data*           modeldata,           /**< model data */
        solutionWindow***     solutionwindows,
        SCIP_Bool*            isfeasible,
        int*                  tour,
        int*                  duration,
        int                   tourlength,
        int                   day
)
{
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    double objective = 0.0;
    int* arrivaltimes;
    int v = modeldata->nC - 1;
    int tourstart;
    int i, j, gamma;
    int tmpduration;
    double windowweight = 0.0;
    modelWindow* window;

    /* start at depot node */
    assert(modeldata->customerIDs[v] == -1);
    assert(modeldata != NULL);
    assert(probdata != NULL);
    assert(tour != NULL);
    assert(tourlength > 0);
    assert(modeldata->maxDelayEvents >= 0);
    if(duration != NULL) *duration = 0;

    SCIP_CALL (SCIPallocMemoryArray(scip, &arrivaltimes, modeldata->maxDelayEvents + 1));
    for (gamma = 0; gamma < modeldata->maxDelayEvents + 1; gamma++)
    {
        arrivaltimes[gamma] = modeldata->shift_start;
    }

    if (solutionwindows != NULL)
    {
        SCIP_CALL( SCIPallocMemoryArray(scip, solutionwindows, tourlength) );
    }
    for (i = 0; i < tourlength; i++)
    {
        assert(v != tour[i]);
        assert(0 <= tour[i] && tour[i] < modeldata->nC - 1);
        /* reuse the memory of arrivaltimes for old and new node --> go backwards through the array and overwrite entries */
        SCIP_CALL( updateRobustTimes(modeldata, arrivaltimes, v, tour[i], day, isfeasible, &windowweight, TRUE, &window));
        if (!(*isfeasible))
        {
            if (solutionwindows != NULL)
            {
                for (j = 0; j < i; j++)
                {
                    SCIP_CALL( freeSolutionWindow(scip, (*solutionwindows)[j]) );
                }
                SCIPfreeMemoryArray(scip, solutionwindows);
            }
            SCIPfreeMemoryArray(scip, &arrivaltimes);
            return 0.0;
        }
        assert(window != NULL);
        if (i == 0)
        {
            tourstart = arrivaltimes[0] - getTravelTime(modeldata, v, tour[i], arrivaltimes[0]);
        }
        /* increase objective value */
        objective += probdata->alphas[1] * getTravelTime(modeldata, v, tour[i], arrivaltimes[0])
                   + probdata->alphas[2] * (1 - windowweight) * WINDOW_WEIGHT_FACTOR;

        if (arrivaltimes[modeldata->maxDelayEvents] - arrivaltimes[0] > probdata->delayTolerance)
        {
            objective += probdata->alphas[0] * modeldata->obj[tour[i]] * (arrivaltimes[modeldata->maxDelayEvents] - arrivaltimes[0] - probdata->delayTolerance);
        }

        //Add price collecting term, if customer was set as optional.
        if (probdata->useOptionals && probdata->optionalCustomers[tour[i]] == TRUE)
        {
            objective -= modeldata->obj[tour[i]] * PRICE_COLLECTING_WEIGHT;
        }

        /* if set, create the solutionWindow for this customer */
        if (solutionwindows != NULL)
        {
            SCIP_CALL( createSolutionWindow(scip, &(*solutionwindows)[i], day, arrivaltimes[0], arrivaltimes[0] + probdata->delayTolerance, arrivaltimes[modeldata->maxDelayEvents], window) );
            assert((*solutionwindows)[i] != NULL);
        }
        v = tour[i];
        assert(v != modeldata->nC -1);
    }
    /* add cost of edge back to depot */
    objective += probdata->alphas[1] * getTravelTime(modeldata, v, modeldata->nC - 1, arrivaltimes[0]);

    SCIP_CALL( updateRobustTimes(modeldata, arrivaltimes, v, modeldata->nC - 1, day, isfeasible, NULL, TRUE, &window));

    if (!(*isfeasible))
    {
        if (solutionwindows != NULL)
        {
            for (j = 0; j < i; j++)
            {
                SCIP_CALL( freeSolutionWindow(scip, (*solutionwindows)[j]) );
            }
            SCIPfreeMemoryArray(scip, solutionwindows);
        }
        SCIPfreeMemoryArray(scip, &arrivaltimes);

        return 0.0;
    }
    if(arrivaltimes[modeldata->maxDelayEvents] - tourstart > WORKTIME_LIMIT)
    {
        *isfeasible = FALSE;
        SCIPfreeMemoryArray(scip, &arrivaltimes);
        return -1;
    }
    tmpduration = arrivaltimes[modeldata->maxDelayEvents] - tourstart;
    if(duration != NULL)
    {
        *duration = tmpduration;
    }
    SCIPfreeMemoryArray(scip, &arrivaltimes);
    return objective;
}

static
int getEdgeArray(
        SCIP*                 scip,
        model_data*           modeldata,
        tuple*                sortededges
){
    int i,j,nC;
    int numedges;
    int** costMatrix;
    int** edgeIndex;
    neighbor* nb;
    int index;

    nC = modeldata->nC;

    SCIP_CALL(SCIPallocMemoryArray(scip, &costMatrix, nC));
    SCIP_CALL(SCIPallocMemoryArray(scip, &edgeIndex, nC));

    for(i = 0; i < nC; i++)
    {
        SCIP_CALL(SCIPallocMemoryArray(scip, &costMatrix[i], nC));
        SCIP_CALL(SCIPallocMemoryArray(scip, &edgeIndex[i], nC));
    }

    for(i = 0; i < nC; i++)
    {
        for(j = i + 1; j < nC; j++)
        {
            edgeIndex[i][j] = 0;
            edgeIndex[j][i] = 0;
            costMatrix[i][j] = -max(4, -modeldata->t_travel[i][j],-modeldata->t_travelAM[i][j],-modeldata->t_travelPM[i][j],-modeldata->t_travelNoon[i][j]);
            costMatrix[j][i] = -max(4, -modeldata->t_travel[j][i],-modeldata->t_travelAM[j][i],-modeldata->t_travelPM[j][i],-modeldata->t_travelNoon[j][i]);
        }
    }
    numedges = 0;
    if(modeldata->neighbors != NULL)
    {
        for(i = 0; i < nC; i++)
        {
            for(j = 0; j < modeldata->nDays; j++)
            {
                nb = modeldata->neighbors[i][j];
                while(nb != NULL)
                {
                    index = edgeIndex[i][nb->id] - 1;
                    if(index == -1)
                    {
                        sortededges[numedges].index = nC * i + nb->id;
                        sortededges[numedges++].value = costMatrix[i][nb->id];
                        edgeIndex[i][nb->id] = numedges;
                        edgeIndex[nb->id][i] = numedges;
                    }else if(sortededges[index].value > costMatrix[i][nb->id])
                    {
                        sortededges[index].value = costMatrix[i][nb->id];
                    }
                    nb = nb->next;
                }
            }
        }
    }else{
        for(i = 0; i < nC; i++)
        {
            for(j = i + 1; j < nC; j++)
            {
                sortededges[numedges].index = nC * i + j;
                sortededges[numedges++].value = -max(2, -costMatrix[i][j], -costMatrix[j][i]);
            }
        }
    }
    for(i = 0; i < nC; i++)
    {
        SCIPfreeMemoryArray(scip, &costMatrix[i]);
        SCIPfreeMemoryArray(scip, &edgeIndex[i]);
    }
    SCIPfreeMemoryArray(scip, &costMatrix);
    SCIPfreeMemoryArray(scip, &edgeIndex);

    return numedges;
}

static
int find(
        int*                  parent,
        int                   u
){
    if(parent[u] != u)
    {
        parent[u] = find(parent, parent[u]);
    }
    return parent[u];
}

static
void merge(
        int*                  parent,
        int                   u,
        int                   v
){
    u = find(parent, u);
    v = find(parent, v);
    parent[v] = u;
}

int getMSTcost(
        SCIP*                 scip,
        model_data*           modeldata
){
    int nC, i, j, e;
    tuple* sortededges;
    int numedges;
    int mstCost;
    int count;
    int* parent;

    nC = modeldata->nC;
    SCIP_CALL(SCIPallocMemoryArray(scip, &sortededges, nC*(nC-1)/2));
    SCIP_CALL(SCIPallocMemoryArray(scip, &parent, nC));

    numedges = getEdgeArray(scip, modeldata, sortededges);

    SCIP_CALL(SCIPreallocMemoryArray(scip, &sortededges, numedges));
    qsort(sortededges, numedges, sizeof(sortededges[0]), cmp_vrp);

    /** Kruskal algorithm */
    mstCost = 0;
    for(i = 0; i < nC; i++)
    {
        parent[i] = i;
    }
    count = 0;
    for(e = numedges - 1; e >= 0; e--)
    {
        i = (int) sortededges[e].index/nC;
        j = (int) sortededges[e].index - (i * nC);
        assert(i != j);
        if(find(parent, i) != find(parent, j))
        {
            count++;
            merge(parent, i, j);
            mstCost += sortededges[e].value;
        }
    }
    assert(count == nC - 1);

    SCIPfreeMemoryArray(scip, &sortededges);
    SCIPfreeMemoryArray(scip, &parent);
    return mstCost;
}

/**
 * Computes a lower bound for the master problem based on simple MST-costs.
 **/
SCIP_RETCODE computeStaticDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   model_data*           modeldata           /**< model data */
   )
{
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   int dualboundMst = 0;
   int dualboundOutgoingArcs = 0;
   double dualbound = 0.0;
   int i, j;
   
   dualboundMst += getMSTcost(scip, modeldata);

   for (i = 0; i < modeldata->nC; i++)
   {
      int lowest = INT_MAX;
      if (modeldata->neighbors != NULL)
      {
         int day;
         for (day = 0; day < modeldata->nDays; day++)
         {
            int lowestday = INT_MAX;
            neighbor* node = modeldata->neighbors[i][day];
            while(node != NULL)
            {
               if (lowest > modeldata->t_travel[i][node->id])
               {
                  lowest = modeldata->t_travel[i][node->id];
               }
               if( lowestday > modeldata->t_travel[i][node->id])
               {
                     lowestday = modeldata->t_travel[i][node->id];
               }
               node = node->next;
            }
            probdata->shortestEdge[day][i] = lowestday;
         }
      } else {
         for (j = 0; j < modeldata->nC; j++)
         {
            if (i == j)
            {
               continue;
            }
            if (lowest > modeldata->t_travel[i][j])
            {
               lowest = modeldata->t_travel[i][j];
            }
         }
      }
      dualboundOutgoingArcs += lowest;
   }

   dualbound = max(2, dualboundMst, dualboundOutgoingArcs);
   dualbound *= probdata->alphas[1];

   /* price collecting */
   if (probdata->useOptionals == TRUE)
   {
      dualbound -= probdata->optionalCost * PRICE_COLLECTING_WEIGHT;
   }
   

   SCIP_CALL( SCIPupdateLocalDualbound(scip, dualbound) );
   return SCIP_OKAY;
}
