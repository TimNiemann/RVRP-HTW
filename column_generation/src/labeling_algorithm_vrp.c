/**@file   labeling_algorithm_vrp.c
 * @brief  labeling Algorithm for column generation 
 * @author Tim Niemann, TU Braunschweig
 */

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <pthread.h>

#include "tools_data.h"
#include "pricer_vrp.h"
#include "probdata_vrp.h"
#include "vardata_vrp.h"
#include "label_vrp.h"
#include "labellist_vrp.h"
#include "tools_vrp.h"
#include "labeling_algorithm_vrp.h"
#include "cons_arcflow.h"

/** returns accumulated lengths of labellists */
static
int labellistsTotalLength(
   label_list**         labellists,
   int                  nlists
   )
{
   int i;
   int length = 0;
   assert(labellists != NULL);
   for (i = 0; i < nlists; i++)
   {
      length += labellistLength(labellists[i]);
   }
   return length;
}

/** Randomly chooses the next customer for propagation */
static
SCIP_RETCODE getNextList(
        SCIP*             scip,
        label_list**      labellists,
        label_list**      list,
        int               nlists
)
{
    int randIndex = rand() % nlists;
    int i;

    assert(labellists != NULL);
    assert(*list == NULL);
    assert(nlists > 0);

    /* search for the next list after randIndex, that is not empty */
    for (i = 0; i < nlists; i++)
    {
        if (labellists[randIndex] != NULL)
        {
            *list = labellists[randIndex];
            labellists[randIndex] = (*list)->next;
            if((*list)->next != NULL)
            {
                (*list)->next->prev = NULL;
            }
            break;
        }
        randIndex = (randIndex + 1) % nlists;
    }
    return SCIP_OKAY;
}

/** Sorts the neighbors of each customer available on this day by dualvalues */
static
SCIP_RETCODE getNeighborsSorted(
   SCIP*                scip,
   model_data*          modeldata,
   tuple**              permutedNeighbors,
   int*                 npermutedNeighbors,
   double*              dualvalues,
   int                  day
   )
{
   int i;
   assert(permutedNeighbors != NULL);
   assert(npermutedNeighbors != NULL);
   assert(dualvalues != NULL);
   for (i = 0; i < modeldata->nC; i++)
   {
      npermutedNeighbors[i] = getNumberOfNeighbors(scip, modeldata, i, day);
      SCIP_CALL( SCIPallocMemoryArray(scip, &(permutedNeighbors[i]), npermutedNeighbors[i]) );

      SCIP_CALL( sortNeighborsOfNode(scip, modeldata, day, i, permutedNeighbors[i], npermutedNeighbors[i], dualvalues));
   }

   return SCIP_OKAY;
}

/** Calculates the latest arrival of each customer of this day */
static
SCIP_RETCODE getUpperTimeWindows(
   SCIP*                scip,
   model_data*          modeldata,
   double*              dualvalues,
   tuple**              permutedNeighbors,
   int*                 npermutedNeighbors,
   int*                 upperTimeWindows,
   int                  day
   )
{
    int i;
    neighbor* nb;
    modelWindow* window;
    assert(dualvalues != NULL);
    assert(permutedNeighbors != NULL);
    assert(npermutedNeighbors != NULL);
    assert(upperTimeWindows != NULL);

    for (i = 0; i < modeldata->nC; i++)
    {
      upperTimeWindows[i] = 0;
    }
    /* just for all customers of this day, the upperTimeWindow will be set to the according nonzero value */
    nb = modeldata->neighbors[modeldata->nC - 1][day];
    while (nb != NULL)
    {
        i = nb->id;
        window = modeldata->timeWindows[i];//getNextTimeWindow(modeldata, i, day, modeldata->shift_start);
        while(window != NULL)
        {
            if(window->day != day)
            {
                window = window->next;
                continue;
            }
            if(window->end_t > upperTimeWindows[i])
            {
                upperTimeWindows[i] = window->end_t;
            }
            if (upperTimeWindows[modeldata->nC - 1] < window->start_t)
            {
                upperTimeWindows[modeldata->nC - 1] = window->start_t;
            }
            window = window->next;
        }
        nb = nb->next;
    }
   return SCIP_OKAY;
}

/** Main method of the labeling algorithm
 * Calculates tours with minimal reduced costs */
static
SCIP_RETCODE generateLabels(
   SCIP* scip,
   model_data* modeldata,
   label_list** bestLabels,
   double* dualvalues,
   SCIP_Bool* visited,
   SCIP_Bool isFarkas,
   SCIP_Bool isHeuristic,
   int day,
   int nUsedNeighbors,
   SCIP_Bool* toDepot
   )
{
   label_list** labellists = NULL;          /* a label list for each customer */
   label_list** usedlabellists = NULL;      /* a label list for each customer for propagated labels */
   label_list* depotlist = NULL;
   labelVrp* label = NULL;
   double sumNegativeRedCosts = 0.0;
   double bestRedCost = 0.0;
   tuple** permutedNeighbors;
   int* npermutedNeighbors;
   int* upperTimeWindows;
   int npropagatedLabels = 0;
   int nbestLabels = 0;
   int nlabels = 0;
   int nUsedLabels = 0;
   int *nUsedLab;
   time_t starttime;
   int i;
    int maxdomi = 0;
    int test;
    int totaldomi = 0;
    int totaldeleted = 0;
    int deletedLabels;
    int dominanceStatus;

   assert(dualvalues != NULL);
   assert(bestLabels != NULL);
   /* if set, the labels with smallest, positive reduced costs are added if there are none with negative cost */
   if (ADD_LABELS_POSITIVE_COST)
   {
      bestRedCost = SCIP_DEFAULT_INFINITY;
   }

   /* Compute the maximum possible reduced costs, every tour could collect */
   sumNegativeRedCosts = - sumOfPossibleDualvalues(scip, modeldata, dualvalues, day, isFarkas);
   assert(sumNegativeRedCosts <= 0);
   if (!SCIPisSumNegative(scip, sumNegativeRedCosts - dualvalues[modeldata->nC - 1 + day]))
   {
      return SCIP_OKAY;
   }
   /* precompute the sequence of neighbors for every node, sorted by dualvalues */
   SCIP_CALL( SCIPallocMemoryArray(scip, &nUsedLab, modeldata->nC - 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &permutedNeighbors, modeldata->nC) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &npermutedNeighbors, modeldata->nC) );

   SCIP_CALL( getNeighborsSorted(scip, modeldata, permutedNeighbors, npermutedNeighbors, dualvalues, day) );
   /* compute the upper limit for a possible arrivaltime at each customer */
   SCIP_CALL( SCIPallocMemoryArray(scip, &upperTimeWindows, modeldata->nC) );
   SCIP_CALL( getUpperTimeWindows(scip, modeldata, dualvalues, permutedNeighbors, npermutedNeighbors, upperTimeWindows, day) );
   /* create the initial, empty label */   
   labelVrpCreateEmpty(scip, &label, modeldata->nC, modeldata->maxDelayEvents + 1, - dualvalues[modeldata->nC - 1 + day], sumNegativeRedCosts, day);
   assert(label != NULL);
   nlabels++;

   /* initialize labellists */
   SCIP_CALL( SCIPallocMemoryArray(scip, &labellists, modeldata->nC - 1) );
   SCIP_CALL( SCIPallocMemoryArray(scip, &usedlabellists, modeldata->nC - 1) );
   for (i = 0; i < modeldata->nC - 1; i++)
   {
      labellists[i] = NULL;
      usedlabellists[i] = NULL;
      nUsedLab[i] = 0;
   }
   label->lhs = dualvalues[modeldata->nC - 1 + day];
   SCIP_CALL( labellistCreate(scip, &depotlist, label, 0));
   labellists[0] = depotlist;
   starttime = time(NULL);

   /** labeling algorithm */
   while (npropagatedLabels < nlabels)
   {
      labelVrp* newLabel = NULL;
      label_list* newList = NULL;
      label_list* lastList = NULL;
      label_list* currentList = NULL;
      int p = 0;  /* Number of neighbors this label was already propagated to */
      assert(nlabels - npropagatedLabels == labellistsTotalLength(labellists, modeldata->nC - 1));
      /* get the next label and propagate it to all neighbors */
      label = NULL;
      SCIP_CALL( getNextList(scip, labellists, &currentList, modeldata->nC - 1) );
      label = currentList->label;
      npropagatedLabels++;
      assert(label != NULL);
      /* propagate this label to all neighbors */
      for (i = 0; i < npermutedNeighbors[label->node]; i++)
      {

         newList = NULL;
         /* if this label was already propagated to many neighbors, skip the other ones */
         if (p >= nUsedNeighbors)
         {
            break;
         }
         /* scip the next customer, if this would create a loop or he was already visited on another day */
         if (permutedNeighbors[label->node][i].index == label->node || (isHeuristic && (visited[permutedNeighbors[label->node][i].index] == TRUE)))
         {
            continue;
         }
         newLabel = NULL;
         labelVrpPropagate(scip, modeldata, label, &newLabel, permutedNeighbors[label->node][i].index, dualvalues[permutedNeighbors[label->node][i].index], isFarkas);
         if (newLabel != NULL)
         {
            /* If in all cases this label would generate no label with better redcost than one already generated, delete it */
            if (!SCIPisSumNegative(scip, newLabel->redcost + newLabel->collactableRedcost - bestRedCost))
            {
               labelVrpFree(scip, &newLabel);
               continue;
            }
            newLabel->collactableRedcost = labelVrpCollactableRedCostTimeDependent(scip, modeldata, newLabel, dualvalues, upperTimeWindows);
            if (!SCIPisSumNegative(scip, newLabel->redcost + newLabel->collactableRedcost - bestRedCost))
            {
               labelVrpFree(scip, &newLabel);
               continue;
            }
            /**** dominance check ****/
            dominanceStatus = labellistDominanceCheck(scip, labellists ,&labellists[newLabel->node], newLabel, FALSE, NULL, NULL);
            deletedLabels = 0;

            test = labellistDominanceCheck(scip, labellists, usedlabellists, newLabel, TRUE, &deletedLabels, nUsedLab);
            nlabels -= deletedLabels;
            totaldeleted += deletedLabels;

            if(dominanceStatus == 0)
            {
                dominanceStatus = test;
            }
            if(dominanceStatus > maxdomi) maxdomi = dominanceStatus;
            switch (dominanceStatus)
            {
            /* this label is dominated, no adding to the pool */
            case -1:
               labelVrpFree(scip, &newLabel);
               totaldomi++;
               break;
            /* dominance check was meaningless, add this label */
            case  0:
               assert(newList == NULL);
               labellistInsertNew(scip, &labellists[newLabel->node], newLabel, &newList, newLabel->redcost);
               assert(newList != NULL);
               p++;
               nlabels++;
               break;
            /* If this label dominates another label already in the list, the other one was deleted, add this label */
            default:
               labellistInsertNew(scip, &labellists[newLabel->node], newLabel, &newList, newLabel->redcost);
               assert(newList != NULL);
               nlabels += 1 - dominanceStatus;
               p++;
               totaldomi += dominanceStatus;
               if(dominanceStatus < -1)
               {
                   SCIPwarningMessage(scip, "Error in dominance check, return value was: %d\n", dominanceStatus);
                   SCIPABORT();
               }
               break;
            }
            if(dominanceStatus >= 0)
            {
                /* if label was added, set label-tree data */
                if(lastList == NULL)
                {
                    currentList->child = newList;
                }else
                {
                    lastList->nextSibling = newList;
                    newList->prevSibling = lastList;
                }
                newList->parent = currentList;
                lastList = newList;
            }
         }
      }

      /* also propagate this label to the depot, if it is not the initial label */
      if (label->node == modeldata->nC -1)
      {
          continue;
      }
      /* continue if the arc to the depot is not available due to branching decisions */
      newLabel = NULL;
      if(toDepot[label->node])
      {
          labelVrpPropagate(scip, modeldata, label, &newLabel, modeldata->nC - 1, 0, isFarkas);
      }

      /* add propagated label to usedlabellist of current node */
      currentList->isPropagated = TRUE;
      currentList->next = NULL;
      currentList->prev = NULL;
      labellistInsertNew(scip, &usedlabellists[label->node], label, &currentList, label->redcost);
      nUsedLab[label->node]++;

      if (newLabel != NULL)
      {
         /* If this is a label with negative reduced costs, which is feasible,
          * add it to the pool and update the value for the current best reduced costs */
         if (SCIPisSumNegative(scip, newLabel->redcost - bestRedCost))
         {
            bestRedCost = newLabel->redcost;
            SCIP_CALL( labellistInsert(scip, bestLabels, newLabel, newLabel->redcost));
            assert(*bestLabels != NULL);
            nbestLabels++;
         } else {
            labelVrpFree(scip, &newLabel);
            newLabel = NULL;
         }
      }

      /* early stopping is not allowed, if this is exact pricing */
      if (!isHeuristic)
      {
         if ((npropagatedLabels % (MAX_CREATED_LABELS * 10)) == 0)
         {
            if (PRINT_EXACT_LABELING)
            {
                nUsedLabels = 0;
                for(i = 0; i < modeldata->nC - 1; i++) nUsedLabels += nUsedLab[i];

               printf("day: %d, nlabels: %d, npropagated: %d, nbest: %d, toCheck: %d, active old Labels: %d, deleted labels: %d\n",
                      day, nlabels, npropagatedLabels, nbestLabels, nlabels-npropagatedLabels, nUsedLabels, totaldeleted+npropagatedLabels-nUsedLabels);
            }
         }
         if (nbestLabels > MAX_CREATED_LABELS || (nlabels * nbestLabels > 500000))
         {
            SCIPdebugMessage("Exact labeling cancelled by too many labels with negative reduced costs.\n");
            if (PRINT_EXACT_LABELING)
            {
               printf("day: %d, nlabels: %d, npropagated: %d, nbest: %d, toCheck: %d\n", day, nlabels, npropagatedLabels, nbestLabels, nlabels-npropagatedLabels);
            }
            break;
         }
         if (time(NULL) - starttime > LABELING_TIME_LIMIT && HEURISTIC_DOMINANCE)
         {
            if (nUsedNeighbors >= modeldata->nC - 1)
            {
               SCIPwarningMessage(scip, "Labeling-Iteration cancelled by time limit on day %d. Optimality of computed solution is not guaranteed.\n", day);
            }
            break;
         }
         continue;
      }

      /* Heuristic call:
      *  pricing is stopped early, if there are too many labels or almost every label has negative reduced cost */
      if (npropagatedLabels > MAX_CREATED_LABELS || (nbestLabels > MIN_REQUIRED_LABELS && (npropagatedLabels/nbestLabels) < 1.5))
      {
         assert(isHeuristic);
         SCIPdebugMessage("Heuristic pricing cancelled by too many labels\n");
         break;
      }
   }

   /* free memory */
   if(depotlist != NULL) labellistDestroy(scip, depotlist);
   for (i = 0; i < modeldata->nC; i++)
   {
      SCIPfreeMemoryArray(scip, &(permutedNeighbors[i]));
   }
   for (i = 0; i < modeldata->nC - 1; i++)
   {
      SCIP_CALL( labellistDestroy(scip, usedlabellists[i]) );
      SCIP_CALL( labellistDestroy(scip, labellists[i]) );
   }
   SCIPfreeMemoryArray(scip, &nUsedLab);
   SCIPfreeMemoryArray(scip, &usedlabellists);
   SCIPfreeMemoryArray(scip, &labellists);
   SCIPfreeMemoryArray(scip, &permutedNeighbors);
   SCIPfreeMemoryArray(scip, &npermutedNeighbors);
   SCIPfreeMemoryArray(scip, &upperTimeWindows);
   
   return SCIP_OKAY;
}

/* add the best label(s) to the master problem */
static
SCIP_RETCODE addToursToMaster(
   SCIP* scip,
   model_data* modeldata,
   label_list** bestLabels,
   SCIP_Bool* visited,
   SCIP_Bool isFarkas,
   int day
   )
   {
   int naddedLabels = 0;
   int nbestLabels = labellistLength(*bestLabels);

   while(naddedLabels < MAX_ADDED_LABELS && nbestLabels > 0)
   {
      SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
      labelVrp* newLabel = NULL;
      char name[SCIP_MAXSTRLEN];
      char strtmp[SCIP_MAXSTRLEN];
      int i;
      assert(bestLabels != NULL);

      SCIP_CALL(labellistExtractFirst(scip, bestLabels, &newLabel));
      nbestLabels--;

      /* double check, that the depot is the last visited node */
      assert(newLabel->visitednodes[newLabel->nvisitednodes - 1] == modeldata->nC -1);

      /* create variable name */
      if (!isFarkas)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingLabelRed_%2d: ", newLabel->day);
      } else {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "pricingLabelFar_%2d: ", newLabel->day);
      }
      for (i = 0; i < newLabel->nvisitednodes - 1; i++)
      {
         (void) SCIPsnprintf(strtmp, SCIP_MAXSTRLEN, "_%d", newLabel->visitednodes[i]);
         strcat(name, strtmp);
      }

      /* Don't add the variable, if it is already contained in prob data*/
      if (!SCIPprobdataContainsVar(probdata, name))
      {
         /* create tour */
         SCIP_Bool isFeasible;
         solutionWindow** solutionwindows = NULL;
         int expectedDuration;
         int tourduration = newLabel->arrivaltimes[newLabel->narrivaltimes - 1] - newLabel->starttime;
         double obj = computeObjValue(scip, modeldata, &solutionwindows, &isFeasible, newLabel->visitednodes, &expectedDuration, newLabel->nvisitednodes - 1, day);
         if(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) == 1)
         {
            /* just for root node, since rearrangeTour does not respect branching decisions */ 
            SCIP_CALL(rearrangeTour(scip, modeldata, newLabel->visitednodes, newLabel->nvisitednodes - 1, &obj, day));
         }
         assert(tourduration == expectedDuration);
         assert(isFeasible);
         assert(solutionwindows != NULL);
         naddedLabels++;

         /* safe visited customers for heuristic call */
         if (visited != NULL)
         {
            for (i = 0; i < newLabel->nvisitednodes - 1; i++)
            {
               visited[newLabel->visitednodes[i]] = TRUE;
            }
         }
         /* Add variable to model */
         SCIP_CALL( SCIPcreateColumn(scip, probdata, name, FALSE, obj, newLabel->visitednodes, newLabel->nvisitednodes - 1, expectedDuration, solutionwindows, day) );
         SCIP_CALL( freeSolutionWindowArray(scip, solutionwindows, newLabel->nvisitednodes - 1) );
      }

      /* free memory */
      labelVrpFree(scip, &newLabel);
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE labelingAlgorithm(
   SCIP*                scip,
   SCIP_Bool            isFarkas,        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
   SCIP_Bool            isHeuristic,
   int                  day,
   SCIP_Bool*           visited,
   SCIP_Bool*           toDepot,
   label_list**         bestLabels
   )
{
   SCIP_PRICER* pricer = NULL;
   SCIP_PRICERDATA* pricerdata = NULL;
   model_data* modeldata = NULL;
   double* dualvalues = NULL;
    int nUsedNeighbors;

   /* get the pricer, problem and model data */
   assert(scip != NULL);
   pricer = SCIPfindPricer(scip, "vrp");
   assert(pricer != NULL);
   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   modeldata = pricerdata->modeldata;
   assert(modeldata != NULL);
   assert(modeldata->nC - 1 + modeldata->nDays == pricerdata->nconss);
   assert((isHeuristic && (visited != NULL)) || (!isHeuristic && (visited == NULL)));

   /* get dual/farkas values */
   SCIP_CALL( SCIPallocMemoryArray(scip, &dualvalues, pricerdata->nconss) );
   SCIP_CALL( getDualValues(scip, dualvalues, isFarkas));

   /* increase the neighborhood size in each iteration */
    nUsedNeighbors = (40 <= modeldata->day_sizes[day] ? 40 : modeldata->day_sizes[day]);
    if (!isHeuristic)
    {
        nUsedNeighbors = (20 <= modeldata->day_sizes[day] ? 20 : modeldata->day_sizes[day]);
    }
   /* generate labels with negative reduced costs and save them in bestLabels */
   while(*bestLabels == NULL && nUsedNeighbors <= modeldata->day_sizes[day])
   {
      nUsedNeighbors*=2;
      SCIP_CALL( generateLabels(scip, modeldata, bestLabels, dualvalues, visited, isFarkas, isHeuristic, day, nUsedNeighbors, toDepot) );
   }
   SCIPfreeMemoryArray(scip, &dualvalues);
   
   return SCIP_OKAY;
}

/** labeling algorithm on one thread */
SCIP_RETCODE labelingAlgorithmIterativ(
        SCIP*                scip,
        SCIP_Bool            isFarkas,        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
        SCIP_Bool            isHeuristic,
        int                  nDays,
        tuple*               days,
        SCIP_Bool*           visited,
        SCIP_Bool*           toDepot
){
    label_list* bestLabels = NULL;
    int i;
    for(i = 0; i < nDays; i++)
    {
        labelingAlgorithm(scip, isFarkas, isHeuristic, days[i].index, visited, toDepot, &bestLabels);

        SCIP_CALL( addToursToMaster(scip, SCIPgetProbData(scip)->modeldata, &bestLabels, visited, isFarkas, days[i].index) );
        SCIP_CALL( labellistDestroy(scip, bestLabels) );
        bestLabels = NULL;
    }

    return SCIP_OKAY;
}

static
void *labeling_thread(void *arguments){
   arg_struct* args = arguments;
   int day = args->day;

   assert(args != NULL);
   assert(!args->isHeuristic && args->visited == NULL);
   
   labelingAlgorithm(args->scip, args->isFarkas, args->isHeuristic, day, args->visited, args->toDepot, &args->bestLabels);

   if (PRINT_EXACT_LABELING)
   {
      printf("Thread for day %d: Ended.\n", day);
   }

   return NULL;
}

/** Same as labeling Algorithm, but uses a different thread for each day */
SCIP_RETCODE labelingAlgorithmParallel(
   SCIP*                scip,
   SCIP_Bool            isFarkas,        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
   SCIP_Bool            isHeuristic,
   int                  nDays,
   SCIP_Bool*           visited,
   SCIP_Bool*           toDepot
   )
{
   pthread_t threads[nDays];
   arg_struct* thread_args;
   int result_code;
   int i;

   SCIP_CALL( SCIPallocMemoryArray(scip, &thread_args, nDays) );

   if (PRINT_EXACT_LABELING)
   {
      printf("Starting Labeling Algorithm in Parallel.\n");
   }
   //create all threads one by one
   for (i = 0; i < nDays; i++) {
      thread_args[i].scip = scip;
      thread_args[i].isFarkas = isFarkas;
      thread_args[i].isHeuristic = isHeuristic;
      thread_args[i].day = i;
      thread_args[i].visited = visited;
      thread_args[i].toDepot = toDepot;
      thread_args[i].bestLabels = NULL;

      result_code = pthread_create(&threads[i], NULL, labeling_thread, &thread_args[i]);
      assert(!result_code);
   }

   //wait for each thread to complete
   for (i = 0; i < nDays; i++) {
      result_code = pthread_join(threads[i], NULL);
      assert(!result_code);
   }

    /* add the best labels as tours of each day to the master problem */
    for(i = 0; i < nDays; i++)
    {
        SCIP_CALL( addToursToMaster(scip, SCIPgetProbData(scip)->modeldata, &(thread_args[i].bestLabels), visited, isFarkas, i) );
        SCIP_CALL( labellistDestroy(scip, thread_args[i].bestLabels) );
    }

   SCIPfreeMemoryArray(scip, &thread_args);

   return SCIP_OKAY;
}