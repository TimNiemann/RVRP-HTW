/**@file   label_vrp.c
 * @brief  label data structure and label propagation function
 * @author Tim Niemann, TU Braunschweig
 */

#include <assert.h>

#include "scip/scip.h"
#include "scip/scipdefplugins.h"

#include "label_vrp.h"
#include "tools_data.h"
#include "tools_vrp.h"
#include "probdata_vrp.h"


/**
 * Local functions
 */

static
SCIP_Bool labelVrpIsFeasible(
    model_data*     modeldata,
    labelVrp*       label
    )
{
    modelWindow* timewindow;
    int i;

    assert(modeldata != NULL);
    assert(label != NULL);

    /* no customer should be visited twice */
    for (i = 0; i < label->nvisitednodes - 1; i++)
    {
        if (label->visitednodes[i] == label->visitednodes[label->nvisitednodes - 1])
        {
            return FALSE;
        }
    }
    
    /* the shift should not be already over */
    if (label->arrivaltimes[label->narrivaltimes - 1] > modeldata->shift_end)
    {
        return FALSE;
    }

    /* maximum worktime must be respected */
    if (label->arrivaltimes[label->narrivaltimes - 1] - label->starttime > WORKTIME_LIMIT)
    {
        return FALSE;
    }

    /* if the node is the depot, there are no time windows */
    if (label->node == modeldata->nC -1)
    {
        return TRUE;
    }
    /* the latest arrival time must be in the same feasible time window as the arrival time*/
    timewindow = getNextTimeWindow(modeldata, label->node, label->day, label->arrivaltimes[0]);
    if (timewindow->start_t > label->arrivaltimes[0] || timewindow->end_t < label->arrivaltimes[label->narrivaltimes - 1])
    {
        return FALSE;
    }

    /* if all conditions are met, this label is feasible */
    return TRUE;
}


/**
 * Interface functions
 */

/** Creates an empty label at the depot node */
extern
SCIP_RETCODE labelVrpCreateEmpty(
    SCIP*           scip, 
    labelVrp**      label,
    int             ncustomers,
    int             narrivaltimes,
    double          initialRedcost,
    double          collactableRedcost,
    int             day
    )
{
    int i;

    assert(label != NULL);
    assert(*label == NULL);
    assert(scip != NULL);
    assert(ncustomers > 1);
    assert(narrivaltimes > 0);

    SCIP_CALL( SCIPallocMemory(scip, label) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &(*label)->arrivaltimes, narrivaltimes));

    (*label)->sizeBitarray = ncustomers/INT_BIT_SIZE + 1;
    SCIP_CALL( SCIPallocMemoryArray(scip, &(*label)->bitVisitednodes, (*label)->sizeBitarray) );

    (*label)->node = ncustomers - 1;
    (*label)->redcost = initialRedcost;
    (*label)->collactableRedcost = collactableRedcost;
    (*label)->visitednodes = NULL;
    (*label)->nvisitednodes = 0;
    (*label)->narrivaltimes = narrivaltimes;
    (*label)->day = day;
    (*label)->starttime = 0;
    (*label)->ndominated = 0;
    (*label)->lhs = 0.0;

    for (i = 0; i < narrivaltimes; i++) 
    {
        (*label)->arrivaltimes[i] = 0;
    }
    for (i = 0; i < (*label)->sizeBitarray; i++)
    {
        (*label)->bitVisitednodes[i] = 0;
    }


    return SCIP_OKAY;
}

/** Creates a label at node */
extern
SCIP_RETCODE labelVrpCreate(
    SCIP*           scip, 
    labelVrp**      label,
    int             node,
    int             ncustomers,
    int*            visitednodes,
    int             nvisitednodes,
    int*            arrivaltimes,
    int             narrivaltimes,
    double          redcost,
    double          collactableRedcost,
    int             starttime,
    int             day
    )
{
    int i;

    assert(scip != NULL);
    assert(label != NULL);
    assert(*label == NULL);
    assert(node >= 0);
    assert(ncustomers > 1);
    assert(arrivaltimes != NULL);
    assert(narrivaltimes > 0);

    SCIP_CALL( labelVrpCreateEmpty(scip, label, ncustomers, narrivaltimes, redcost, collactableRedcost, day) );

    SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*label)->visitednodes, visitednodes, nvisitednodes) );

    (*label)->node = node;
    (*label)->nvisitednodes = nvisitednodes;
    (*label)->narrivaltimes = narrivaltimes;
    (*label)->day = day;
    (*label)->starttime = starttime;

    for (i = 0; i < narrivaltimes; i++) 
    {
        (*label)->arrivaltimes[i] = arrivaltimes[i];
    }

    /* set bit array of visitednodes */
    for (i = 0; i < nvisitednodes; i++)
    {
        SetBit((*label)->bitVisitednodes, (*label)->visitednodes[i]);
    }

    return SCIP_OKAY;
}

/** Free label data */
extern
SCIP_RETCODE labelVrpFree(
    SCIP*           scip, 
    labelVrp**      label
    )
{
    assert(scip != NULL);
    assert(label != NULL);
    assert(*label != NULL);

    if ((*label)->visitednodes != NULL)
    {
        SCIPfreeMemoryArray(scip, &(*label)->visitednodes);
    }
    SCIPfreeMemoryArray(scip, &(*label)->arrivaltimes);
    SCIPfreeMemoryArray(scip, &(*label)->bitVisitednodes);
    SCIPfreeMemory(scip, label);

    return SCIP_OKAY;
}

/** Propagate a label to destination */
extern
SCIP_RETCODE labelVrpPropagate(
    SCIP*           scip,
    model_data*     modeldata,
    labelVrp*       oldLabel,
    labelVrp**      newLabel,
    int             end,
    double          dualvalue,
    SCIP_Bool       isFarkas
    )
{
    /* potential new label variables */
    int* visitednodes;
    int nvisitednodes;
    int* arrivaltimes;
    int narrivaltimes;
    double redcost;
    double collactableRedcost;
    int starttime;
    int day;
    int start;
    int i;
    SCIP_Bool isfeasible;
    double windowweight = 0.0;
    double lhs;
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

    /* checklist (only in debug mode) */
    assert(scip != NULL);
    assert(modeldata != NULL);
    assert(oldLabel != NULL);
    assert(newLabel != NULL);
    assert(*newLabel == NULL);
    assert(0 <= end && end < modeldata->nC);
    assert(oldLabel->node != end);
    assert(probdata != NULL);

    /* new visited nodes */
    nvisitednodes = oldLabel->nvisitednodes + 1;
    SCIP_CALL( SCIPallocMemoryArray(scip, &visitednodes, nvisitednodes));
    for (i = 0; i < oldLabel->nvisitednodes; i++) 
    {
        visitednodes[i] = oldLabel->visitednodes[i];
    }
    visitednodes[nvisitednodes - 1] = end;

    start = oldLabel->node;
    day = oldLabel->day;
    narrivaltimes = oldLabel->narrivaltimes;
    SCIP_CALL( SCIPduplicateMemoryArray(scip, &arrivaltimes, oldLabel->arrivaltimes,narrivaltimes) );

    /* Calculate new arrival times */
    SCIP_CALL( updateRobustTimes(modeldata, arrivaltimes, start, end, day, &isfeasible, &windowweight, FALSE, NULL));

    if (!isfeasible && end != modeldata->nC - 1)
    {
        SCIPfreeMemoryArray(scip, &visitednodes);
        SCIPfreeMemoryArray(scip, &arrivaltimes);
        return SCIP_OKAY;
    }

    /* (new) start time */
    if (nodeIsDepot(modeldata, oldLabel->node))
    {
        starttime = arrivaltimes[0] - getTravelTime(modeldata, modeldata->nC - 1, end, arrivaltimes[0]);
        assert(starttime >= modeldata->shift_start);
    } else {
        starttime = oldLabel->starttime;
    }

    /* new reduced costs */
    redcost = oldLabel->redcost - dualvalue;
    lhs = oldLabel->lhs + dualvalue;
    collactableRedcost = oldLabel->collactableRedcost;

    /* if this is no farkas pricing, add the objective function to the red cost */
    if (!isFarkas)
    {
        assert(probdata != NULL);
        redcost += probdata->alphas[1] * getTravelTime(modeldata, start, end, arrivaltimes[0])
                +  probdata->alphas[2] * windowweight * WINDOW_WEIGHT_FACTOR;
        /* delay objective */
        if (arrivaltimes[narrivaltimes - 1] - arrivaltimes[0] > probdata->delayTolerance)
        {
            redcost += probdata->alphas[0] * modeldata->obj[end] * (arrivaltimes[narrivaltimes - 1] - arrivaltimes[0] - probdata->delayTolerance);
        }
        /* Price Collecting for hard customers */
        if (probdata->useOptionals == TRUE && probdata->optionalCustomers[end] == TRUE)
        {
            redcost -= modeldata->obj[i] * PRICE_COLLECTING_WEIGHT;
        }
    }

    /* check if there are enough collectable dual values left to create a tour with negative reduced costs */
    if (SCIPisSumNegative(scip, redcost + collactableRedcost) && !SCIPisSumPositive(scip, collactableRedcost))
    {
        /* create the new label */
        SCIP_CALL( labelVrpCreate(scip, newLabel, end, modeldata->nC, visitednodes, nvisitednodes, arrivaltimes, narrivaltimes, redcost, collactableRedcost, starttime, day) );

        (*newLabel)->lhs = lhs;

        for (i = 0; i < (*newLabel)->sizeBitarray; i++)
        {
            (*newLabel)->bitVisitednodes[i] = oldLabel->bitVisitednodes[i];
        }
        SetBit((*newLabel)->bitVisitednodes, end);
        /* check if the newly created label is feasible */
        if (!labelVrpIsFeasible(modeldata, *newLabel))
        {
            labelVrpFree(scip, newLabel);
            *newLabel = NULL;
        }
    }

    /* free used memory */
    SCIPfreeMemoryArray(scip, &visitednodes);
    SCIPfreeMemoryArray(scip, &arrivaltimes);

    return SCIP_OKAY;
}

/**
 * Checks if the first label dominates the second one
 * @param labelA first label
 * @param labelB second label
 * @return TRUE if labelA dominates labelB, 
 *         FALSE else */
SCIP_Bool labelVrpDominates(
    SCIP*       scip,
    labelVrp*   labelA,
    labelVrp*   labelB
    )
{
    int i;
    assert(labelA != NULL);
    assert(labelB != NULL);

    /* check reduced costs */
    if (SCIPisSumPositive(scip, labelA->redcost - labelB->redcost))
    {
        return FALSE;
    }

    /* check number of visited customers */
    if (labelA->nvisitednodes > labelB->nvisitednodes)
    {
        return FALSE;
    }

    /* start time of labelA must be later */
    if (labelA->starttime < labelB->starttime)
    {
        return FALSE;
    }

    /* the arrival times of labelA must be earlier */
    for (i = 0; i < labelA->narrivaltimes; i++)
    {
        if (labelA->arrivaltimes[i] > labelB->arrivaltimes[i])
        {
            return FALSE;
        }
    }

    /* if set, this is a heuristic dominance check without subset condition */
    if (HEURISTIC_DOMINANCE)
    {
        return TRUE;
    }

    /* the visited nodes of labelA must be a subset of labelB */
    for (i = 0; i < labelA->sizeBitarray; i++)
    {
        if(!IsSubsetBit(labelA->bitVisitednodes[i], labelB->bitVisitednodes[i]))
        {
            return FALSE;
        }
    }

   /* the difference between two successive arrival times must not be higher */
   for(i = 0; i < labelA->narrivaltimes - 1; i++)
   {
       if (labelA->arrivaltimes[i + 1] - labelA->arrivaltimes[i] >
           labelB->arrivaltimes[i + 1] - labelB->arrivaltimes[i])
       {
           return FALSE;
       }
   }

    /* the difference between two successive arrival times must not be higher */
    for (i = 0; i < labelA->narrivaltimes - 1; i++)
    {
        if (labelA->arrivaltimes[i + 1] - labelA->arrivaltimes[i] >
            labelB->arrivaltimes[i + 1] - labelB->arrivaltimes[i])
        {
            return FALSE;
        }
    }

    return TRUE;
}

/**
 * computes still collectable dualvalues for a given label
 * @param label label with arivaltimes and current day
 * @param dualvalues dualvalues
 * @param upperTimeWindows nC-array with the last possible arrivaltime for each customer, 
 *                         this equals upper boundary of time window when customer has time window on this day
 *                         in all other cases the value is zero
 * @return negative sum over all reduced cost that still could be collected for this label*/
double labelVrpCollactableRedCostTimeDependent(
    SCIP*           scip,
    model_data*     modeldata,
    labelVrp*       label,
    double*         dualvalues,
    int*            upperTimeWindows
    )
{
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    double possibleDualvalue = 0.0;
    assert(scip != NULL);
    assert(modeldata != NULL);
    assert(label != NULL);
    assert(dualvalues != NULL);
    assert(upperTimeWindows != NULL);

    if (modeldata->neighbors != NULL)
    {
        neighbor* node = modeldata->neighbors[modeldata->nC - 1][label->day];
        assert(node != NULL);
        while (node != NULL)
        {
            if (node->id == modeldata->nC - 1)
            {
                node = node->next;
                continue;
            }
            if (!TestBit(label->bitVisitednodes, node->id))
            {
                if (label->arrivaltimes[label->narrivaltimes - 1] + getTravelTime(modeldata, label->node, node->id, label->arrivaltimes[label->narrivaltimes - 1]) + modeldata->t_service[label->node] < upperTimeWindows[node->id])
                {
                    double serviceThreshold = probdata->alphas[1] * probdata->shortestEdge[label->day][node->id];
                    /* Price Collecting for hard customers */
                    double hardCustomerBonus = 0;
                    if (probdata->useOptionals == TRUE && probdata->optionalCustomers[node->id] == TRUE)
                    {
                        hardCustomerBonus = modeldata->obj[node->id] * PRICE_COLLECTING_WEIGHT;
                    }
                    /* if the current time is after the beginning of all time windows, there can't be any more resets of delay, so the current delay will be also present in all future nodes */
                    if (label->arrivaltimes[0] >= upperTimeWindows[modeldata->nC - 1] || HEURISTIC_COLLECTABLE)
                    {
                        double delayThreshold = probdata->alphas[0] * modeldata->obj[node->id] * (label->arrivaltimes[label->narrivaltimes - 1] - label->arrivaltimes[0] - probdata->delayTolerance);
                        if (SCIPisSumPositive(scip, dualvalues[node->id] + hardCustomerBonus - serviceThreshold - delayThreshold))
                        {
                            possibleDualvalue -= dualvalues[node->id] + hardCustomerBonus - serviceThreshold - delayThreshold;
                        }
                    } else {
                        if (SCIPisSumPositive(scip, dualvalues[node->id] + hardCustomerBonus - serviceThreshold))
                        {
                            possibleDualvalue -= dualvalues[node->id] + hardCustomerBonus - serviceThreshold;
                        }
                    }
                } else {
                    SetBit(label->bitVisitednodes, node->id);
                }
            }
            node = node->next;
        }
    /* if the instance is not preprocessed */
    } else {
        int i;
        for (i = 0; i < modeldata->nC - 1; i++)
        {
            if (!TestBit(label->bitVisitednodes, i)
                && SCIPisSumPositive(scip, dualvalues[i] - probdata->alphas[1] * probdata->shortestEdge[label->day][i])
                && label->arrivaltimes[label->narrivaltimes - 1] < upperTimeWindows[i])
            {
                possibleDualvalue -= dualvalues[i] - probdata->alphas[1] * probdata->shortestEdge[label->day][i];
            }
        }
    }

    return possibleDualvalue;
}

/** Print label data */
extern
void labelVrpPrint(
    labelVrp*      label
    )
{
    int i;
    assert(label != NULL);
    printf("Node: %d\n", label->node);
    printf("RedCost: %f\n", label->redcost);
    printf("Collactable: %f\n", label->collactableRedcost);
    printf("Starttime: %d\n", label->starttime);
    printf("Duration: %d\n", label->arrivaltimes[label->narrivaltimes - 1] - label->starttime);
    printf("Day: %d\n", label->day);
    printf("#visitednodes: %d\n", label->nvisitednodes);
    printf("#arrivaltimes: %d\n", label->narrivaltimes);
    assert(label->node >= 0);
    assert(label->day >= 0);
    assert(label->nvisitednodes >= 0);
    assert(label->narrivaltimes >= 0);

    for (i = 0; i < label->nvisitednodes; i++) {
        printf("%d-> ", label->visitednodes[i]);
        assert(0 <= label->visitednodes[i] );
    }
    printf("\n");

    for (i = 0; i < label->narrivaltimes; i++) {
        printf("%d, ", label->arrivaltimes[i]);
        assert(0 <= label->arrivaltimes[i] );
    }
    printf("\n\n");
}