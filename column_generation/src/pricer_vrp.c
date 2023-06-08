/**@file   pricer_vrp.c
 * @brief  pricer implementation with an exact labeling algorithm
 * @author Tim Niemann, TU Braunschweig
 */

#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <time.h>

#include "scip/cons_setppc.h"
#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

#include "tools_data.h"
#include "probdata_vrp.h"
#include "pricer_vrp.h"
#include "tools_vrp.h"
#include "labeling_algorithm_vrp.h"
#include "pricing_heuristic_vrp.h"

#include "arcflow_branching.h"
#include "cons_arcflow.h"
#include "primal_heuristic_vrp.h"

/**@name Pricer properties
 *
 * @{
 */

#define PRICER_NAME            "vrp"
#define PRICER_DESC            "pricer for vrp tours"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */

/**@} */

/** Clears the neighborhood at the current node */
static
SCIP_RETCODE clearNeighbors(
    SCIP*              scip,
    SCIP_PRICERDATA*   pricerdata
)
{
    int i, j;
    neighbor* nb;

    assert(pricerdata != NULL);
    assert(pricerdata->neighbors != NULL);

    for(i = 0; i < pricerdata->nC; i++)
    {
        assert(pricerdata->neighbors[i] != NULL);

        for(j = 0; j < pricerdata->nDays; j++)
        {
            while(pricerdata->neighbors[i][j] != NULL)
            {
                nb = pricerdata->neighbors[i][j];
                pricerdata->neighbors[i][j] = pricerdata->neighbors[i][j]->next;
                SCIPfreeBlockMemory(scip, &nb);
            }
        }
    }
    return SCIP_OKAY;
}

/** Calculates the graph of the current branching node
 *  Includes local branching decisions */
static
SCIP_RETCODE getCurrentNeighborhood(
    SCIP*               scip,
    SCIP_PRICERDATA*    pricerdata
)
{
    model_data* modeldata;
    SCIP_CONSHDLR* conshdlr;
    SCIP_CONS** conss;
    SCIP_CONS* cons;
    int ncons;
    int tail, head;
    CONSTYPE type;
    neighbor* nb;
    neighbor* newnb;
    neighbor* lastnb;

    int* successor;
    int* predecessor;
    int i, j;

    assert( pricerdata != NULL );

    modeldata = pricerdata->modeldata;
    conshdlr = pricerdata->conshdlr;

    assert( modeldata != NULL );
    assert( scip != NULL );
    assert( conshdlr != NULL );

    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &successor, modeldata->nC) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &predecessor, modeldata->nC) );
    SCIP_CALL( SCIPallocBlockMemory(scip, &newnb) );

    /* clear old neighborhood */
    if(pricerdata->lastnodeid != 1)
    {
        clearNeighbors(scip, pricerdata);
    }
    else
    {
        for(i = 0; i < pricerdata->modeldata->nC; i++)
        {
            for(j = 0; j < pricerdata->modeldata->nDays; j++)
            {
                pricerdata->neighbors[i][j] = NULL;
            }
        }
    }

    /* safe prohibited and enforced arcs */
    for(i = 0; i < modeldata->nC; i++)
    {
        successor[i] = -1;
        predecessor[i] = -1;
        for(j = 0; j < modeldata->nC; j++)
        {
            pricerdata->isForbidden[i][j] = FALSE;
        }
    }

    /* collect all branching decision constraints */
    conss = SCIPconshdlrGetConss(conshdlr);
    ncons = SCIPconshdlrGetNConss(conshdlr);
    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];

        if(!SCIPconsIsActive(cons)) continue;

        /* branching decision data */
        tail = SCIPgetTailArcflow(scip, cons);
        head = SCIPgetHeadArcflow(scip, cons);
        type = SCIPgetTypeArcflow(scip, cons);

        if( type == PROHIBIT )
        {
            /* arc should not already be forbidden */
            assert(!pricerdata->isForbidden[tail][head]);

            pricerdata->isForbidden[tail][head] = TRUE;
            if(head == modeldata->nC - 1)
            {
                pricerdata->toDepot[tail] = FALSE;
            }
        }
        else if( type == ENFORCE )
        {
            /* each customer can only have one in- and one outgoing enforced arc */
            assert(tail == modeldata->nC - 1 || successor[tail] == -1);
            assert(head == modeldata->nC - 1 || predecessor[head] == -1);

            /* save for neighborhood generation */
            successor[tail] = head;
            predecessor[head] = tail;

            /* calculate induced forbidden arcs */
            if(head != modeldata->nC - 1 && tail != modeldata->nC - 1)
            {
                pricerdata->toDepot[tail] = FALSE;
                /* Safe forbidden arcs */
                for(j = 0; j < modeldata->nC; j++)
                {
                    if(j != head) pricerdata->isForbidden[tail][j] = TRUE;
                    if(j != tail) pricerdata->isForbidden[j][head] = TRUE;
                }
            }else if(head == modeldata->nC - 1)
            {
                for(j = 0; j < modeldata->nC - 1; j++)
                {
                    pricerdata->isForbidden[tail][j] = TRUE;
                }
            }else
            {
                assert(tail == modeldata->nC - 1);
                for(j = 0; j < modeldata->nC - 1; j++)
                {
                    pricerdata->isForbidden[j][head] = TRUE;
                }
            }
        }
        else
        {
            SCIPerrorMessage("unknow constraint type <%d>\n, type");
            return SCIP_INVALIDDATA;
        }
    }

    /* create current neighborhood */
    for(i = 0; i < modeldata->nC - 1; i++)
    {
        for(j = 0; j < modeldata->nDays; j++)
        {
            assert(pricerdata->neighbors[i][j] == NULL);

            /* if customer i has an enforced outgoing arc only set one neighbor */
            if(successor[i] != -1)
            {
                newnb->id = successor[i];
                newnb->next = NULL;
                pricerdata->neighbors[i][j] = newnb;
                SCIP_CALL( SCIPallocBlockMemory(scip, &newnb) );
                continue;
            }
            else
            {
                lastnb = NULL;
                nb = modeldata->neighbors[i][j];
                /* copy from modeldata but skip forbidden neighbors */
                while(nb != NULL)
                {
                    if(pricerdata->isForbidden[i][nb->id] || predecessor[nb->id] != -1)
                    {
                        assert(predecessor[nb->id] != i);
                    }
                    else
                    {
                        newnb->id = nb->id;
                        newnb->next = NULL;
                        if(lastnb == NULL)
                        {
                            pricerdata->neighbors[i][j] = newnb;
                            lastnb = newnb;
                        }
                        else
                        {
                            lastnb->next = newnb;
                            lastnb = newnb;
                        }
                        SCIP_CALL( SCIPallocBlockMemory(scip, &newnb) );
                    }
                    nb = nb->next;
                }
            }
        }
    }
    /* special case - depot */
    i = modeldata->nC - 1;
    for(j = 0; j < modeldata->nDays; j++)
    {
        assert(pricerdata->neighbors[i][j] == NULL);
        lastnb = NULL;
        nb = modeldata->neighbors[i][j];
        while(nb != NULL)
        {
            if(pricerdata->isForbidden[i][nb->id] || (predecessor[nb->id] != i && predecessor[nb->id] != -1))
            {
                nb = nb->next;
                continue;
            }
            else
            {
                newnb->id = nb->id;
                newnb->next = NULL;
                if(lastnb == NULL)
                {
                    pricerdata->neighbors[i][j] = newnb;
                    lastnb = newnb;
                }
                else
                {
                    lastnb->next = newnb;
                    lastnb = newnb;
                }
                SCIP_CALL( SCIPallocBlockMemory(scip, &newnb) );
            }
            nb = nb->next;
        }

    }

    SCIPfreeBlockMemory(scip, &newnb);

    SCIPfreeBlockMemoryArray(scip, &predecessor, modeldata->nC);
    SCIPfreeBlockMemoryArray(scip, &successor, modeldata->nC);

    return SCIP_OKAY;
}


/**name Callback methods
 *
 * @{
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreeVrp)
{
   SCIP_PRICERDATA* pricerdata;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);

   if( pricerdata != NULL)
   {
      /* free memory */
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->conss, pricerdata->nconss);
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->toDepot, pricerdata->nC - 1);
       if(pricerdata->neighbors != NULL)
       {
           if(pricerdata->lastnodeid > 1)
           {
               clearNeighbors(scip, pricerdata);
           }
           for(int i = 0; i < pricerdata->nC; i++)
           {
               SCIPfreeBlockMemoryArray(scip, &pricerdata->neighbors[i], pricerdata->nDays);
           }
           SCIPfreeBlockMemoryArray(scip, &pricerdata->neighbors, pricerdata->nC);
       }
       for(int i = 0; i < pricerdata->nC; i++)
       {
           SCIPfreeBlockMemoryArray(scip, &pricerdata->isForbidden[i], pricerdata->nC);
       }
       SCIPfreeBlockMemoryArray(scip, &pricerdata->isForbidden, pricerdata->nC);

      SCIPfreeBlockMemory(scip, &pricerdata);
   }
   //printf("Deinitialization of pricer successful!\n");
   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitVrp)
{  
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS* cons;
   int c;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get transformed constraints */
   for( c = 0; c < pricerdata->nconss; ++c )
   {
      cons = pricerdata->conss[c];

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->conss[c]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->conss[c]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->conss[c]) );
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolVrp)
{
   SCIP_PRICERDATA* pricerdata;
   int c;

   assert(scip != NULL);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* get release constraints */
   for( c = 0; c < pricerdata->nconss; ++c )
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->conss[c])) );
   }

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostVrp)
{
   tuple* days;
   SCIP_PRICERDATA* pricerdata = SCIPpricerGetData(pricer);
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   SCIP_Bool* visited;
   int i;
   int nvars;

   assert(pricerdata != NULL);
   *result = SCIP_SUCCESS;
   if (SCIPgetDualbound(scip) == -SCIP_DEFAULT_INFINITY)
   {
      SCIP_CALL( computeStaticDualbound(scip, probdata->modeldata));
   }
   /* if set, there is no reduced cost pricing, just farkas pricing to find a feasible solution */
   if (NO_REDCOST_PRICING)
   {
      return SCIP_OKAY;
   }

   if (STOP_IF_FEASIBLE && SCIPgetNSols(scip) > 0)
   {
      return SCIP_OKAY;
   }

    /* If this is the first iteration at the current branching node, set neighborhood */
    if(pricerdata->lastnodeid != SCIPnodeGetNumber(SCIPgetCurrentNode(scip)))
    {
        SCIP_CALL( getCurrentNeighborhood(scip, pricerdata) );

        pricerdata->lastnodeid = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
    }

    /* run primal heuristic, if we found a better LP-solution since last call */
    if(SCIPgetNNodes(scip) >= 1)
    {
        if(SCIPisSumNegative(scip, SCIPgetLPObjval(scip) - SCIPgetSolOrigObj(scip, SCIPgetBestSol(scip)) )
           && (!SCIPisSumEQ(scip, SCIPgetLPObjval(scip), pricerdata->lastLPVal)))
        {
            SCIP_CALL(dispatchingHeuristic(scip, probdata, pricerdata));
            pricerdata->lastLPVal = SCIPgetLPObjval(scip);
        }
    }
    /* Try local search pricing at root node */
    if(SCIPgetNNodes(scip) == 1)
    {
       nvars = SCIPgetNVars(scip);
       heuristicPricing(scip, FALSE);
       if(nvars < SCIPgetNVars(scip))
       {
           return SCIP_OKAY;
       }
   }

   /* In heuristic calls we do not search for multiple tours that visit the same customer */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &visited, pricerdata->modeldata->nC - 1));
   for (i = 0; i < pricerdata->modeldata->nC - 1; i++)
   {
      visited[i] = FALSE;
   }
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays));
   /* Create tuples of days and their respective dualvalue */
   for (i = 0; i < pricerdata->modeldata->nDays; i++)
   {
      days[i].index = i;
      days[i].value = SCIPgetDualsolSetppc(scip, pricerdata->conss[pricerdata->modeldata->nC - 1 + i]);
   }

   /* Sort the tuples by dualvalue, days with large dualvalue come first because this result in smaller reduced costs */
   qsort(days, pricerdata->modeldata->nDays, sizeof(days[0]), cmp_vrp);

   /* first try to find a tour with some heuristics */
   SCIP_CALL( labelingAlgorithmIterativ(scip, FALSE, TRUE, pricerdata->modeldata->nDays, days, visited, pricerdata->toDepot) );

   SCIPfreeBlockMemoryArray(scip, &visited, pricerdata->modeldata->nC - 1);
   /* Success? */
   if(nvars < SCIPgetNVars(scip))
   {
       SCIPfreeBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays);
       return SCIP_OKAY;
   }

   visited = NULL;

   /* if no tour was generated, try exact labeling */

   if (PRINT_EXACT_LABELING)
   {
      printf("Reduced cost pricing: Heuristic unsuccessful, trying exact pricing now.\n");
   }
   if (PARALLEL_LABELING)
   {
      SCIP_CALL( labelingAlgorithmParallel(scip, FALSE, FALSE, pricerdata->modeldata->nDays, NULL, pricerdata->toDepot) );
   } else {
      SCIP_CALL( labelingAlgorithmIterativ(scip, FALSE, FALSE, pricerdata->modeldata->nDays, days, NULL, pricerdata->toDepot) );
   }


   SCIPfreeBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays);
   return SCIP_OKAY;
}


/** farkas pricing method of variable pricer for infeasible LPs */
static
SCIP_DECL_PRICERFARKAS(pricerFarkasVrp)
{
   int i, nvars;
   tuple* days;
   SCIP_PRICERDATA* pricerdata;
   SCIP_Bool* visited;

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);
   *result = SCIP_SUCCESS;
   if (SCIPgetDualbound(scip) == -SCIP_DEFAULT_INFINITY)
   {
      SCIP_CALL( computeStaticDualbound(scip, pricerdata->modeldata));
   }
   if (STOP_IF_FEASIBLE && SCIPgetNSols(scip) > 0)
   {
      return SCIP_OKAY;
   }

   /* Try local search pricing at root node */
   if(SCIPgetNNodes(scip) == 1)
   {
       nvars = SCIPgetNVars(scip);

       heuristicPricing(scip, TRUE);
       if(nvars < SCIPgetNVars(scip))
       {
           return SCIP_OKAY;
       }
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &visited, pricerdata->modeldata->nC - 1));
   for (i = 0; i < pricerdata->modeldata->nC - 1; i++)
   {
      visited[i] = FALSE;
   }
   nvars = SCIPgetNVars(scip);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays));
   /* Create tuples of days and their respective dualvalue */
   for (i = 0; i < pricerdata->modeldata->nDays; i++)
   {
      days[i].index = i;
      days[i].value = SCIPgetDualfarkasSetppc(scip, pricerdata->conss[pricerdata->modeldata->nC - 1 + i]);
   }

   /* Sort the tuples by dualvalue, days with large dualvalue come first because this result in smaller reduced costs */
   qsort(days, pricerdata->modeldata->nDays, sizeof(days[0]), cmp_vrp);

   /* If this is the first iteration at the current branching node, set neighborhood */
   if(pricerdata->lastnodeid != SCIPnodeGetNumber(SCIPgetCurrentNode(scip)))
   {
       SCIP_CALL( getCurrentNeighborhood(scip, pricerdata) );

       pricerdata->lastnodeid = SCIPnodeGetNumber(SCIPgetCurrentNode(scip));
   }

   /* first try to find a tour with some heuristics */
   SCIP_CALL( labelingAlgorithmIterativ(scip, TRUE, TRUE, pricerdata->modeldata->nDays, days, visited, pricerdata->toDepot) );

   SCIPfreeBlockMemoryArray(scip, &visited, pricerdata->modeldata->nC - 1);
   /* success? */
   if(nvars < SCIPgetNVars(scip))
   {
       SCIPfreeBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays);
       return SCIP_OKAY;
   }

   /* if no tour was generated, try exact labeling */

   if (PRINT_EXACT_LABELING)
   {
      printf("Farkas pricing: Heuristic unsuccessful, trying exact pricing now.\n");
   }
   for (i = 0; i < pricerdata->modeldata->nDays; i++)
   {
      nvars = SCIPgetNVars(scip);
      SCIP_CALL( labelingAlgorithmIterativ(scip, TRUE, FALSE, pricerdata->modeldata->nDays, days, NULL, pricerdata->toDepot) );
   }

   SCIPfreeBlockMemoryArray(scip, &days, pricerdata->modeldata->nDays);
   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates the vrp variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerVrp(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create vrp variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );

   pricerdata->conss = NULL;
   pricerdata->nconss = 0;
   pricerdata->conshdlr = SCIPfindConshdlr(scip, "arcflow");

   assert(pricerdata->conshdlr != NULL);
   pricerdata->lastnodeid = 1;
   pricerdata->neighbors = NULL;
   pricerdata->toDepot = NULL;
   pricerdata->nC = 0;
   pricerdata->nDays = 0;
   pricerdata->lastLPVal = DBL_MAX;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostVrp, pricerFarkasVrp, pricerdata) );

   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreeVrp) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitVrp) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolVrp) );


   return SCIP_OKAY;
}

/** added problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerVrpActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< set partitioning constraints for the customers */
   int                   nconss,             /**< number of constraints */
   model_data*           modeldata           /**< model data */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int i, c;

   assert(scip != NULL);
   assert(conss != NULL);
   assert(modeldata != NULL);
   assert(nconss > 0);

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert(pricer != NULL);

   pricerdata = SCIPpricerGetData(pricer);
   assert(pricerdata != NULL);

   /* copy arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->conss, conss, nconss) );
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->neighbors, modeldata->neighbors, modeldata->nC) );
   for(i = 0; i < modeldata->nC; i++)
   {
       SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(pricerdata->neighbors[i]), modeldata->neighbors[i], modeldata->nDays) );
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &pricerdata->toDepot, modeldata->nC - 1) );
   for(i = 0; i < modeldata->nC - 1; i++)
   {
       pricerdata->toDepot[i] = TRUE;
   }
   pricerdata->nC = modeldata->nC;
   pricerdata->nDays = modeldata->nDays;
   pricerdata->nconss = nconss;
   pricerdata->modeldata = modeldata;

   /* capture all constraints */
   for( c = 0; c < nconss; ++c )
   {
      SCIP_CALL( SCIPcaptureCons(scip, conss[c]) );
   }
    /* Allocate memory for forbidden arcs */
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pricerdata->isForbidden), pricerdata->nC));
    for(i = 0; i < pricerdata->nC; i++)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(pricerdata->isForbidden[i]), pricerdata->nC) );
        for(c = 0; c < pricerdata->nC; c++)
        {
            pricerdata->isForbidden[i][c] = FALSE;
        }
    }

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}

/** Sets the current dual values */
extern
SCIP_RETCODE getDualValues(
   SCIP* scip,
   double* dualvalues,
   SCIP_Bool isFarkas
   )
{
   SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
   model_data* modeldata = probdata->modeldata;
   SCIP_CONS** conss = probdata->conss;
   int i, day;
   assert(dualvalues != NULL);

   for (i = 0; i < modeldata->nC - 1; i++)
   {
      assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(conss[i]) ), "setppc", 6) );
      if (isFarkas)
      {
         /* use the farkas value */
         dualvalues[i] = SCIPgetDualfarkasSetppc(scip, conss[i]);
         SCIPdebugMsg(scip, "farkasvalue %d: %f\n", i, dualvalues[i]);
      } else {
         /* or get dual solution vector */
         dualvalues[i] = SCIPgetDualsolSetppc(scip, conss[i]);
         SCIPdebugMsg(scip, "dualvalue %d: %f\n", i, dualvalues[i]);
      }
   }
   for (day = 0; day < modeldata->nDays; day++)
   {
        assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(conss[modeldata->nC - 1 + day]) ), "setppc", 6) );
        if (isFarkas)
        {
            /* use the farkas value */
            dualvalues[modeldata->nC - 1 + day] = SCIPgetDualfarkasSetppc(scip, conss[modeldata->nC - 1 + day]);
            SCIPdebugMsg(scip, "farkasvalue %d: %f\n", i, dualvalues[modeldata->nC - 1 + day]);
        } else {
            /* or get dual solution vector */
            dualvalues[modeldata->nC - 1 + day] = SCIPgetDualsolSetppc(scip, conss[modeldata->nC - 1 + day]);
            SCIPdebugMsg(scip, "dualvalue %d: %f\n", i, dualvalues[modeldata->nC - 1 + day]);
        }
    }

   return SCIP_OKAY;
}

/**
 * @return minimum reduced cost that could be collected */
extern
double sumOfPossibleDualvalues(
    SCIP*           scip,
    model_data*     modeldata,
    double*         dualvalues,
    int             day,
    SCIP_Bool       isFarkas
    )
{
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    modelWindow* window;
    neighbor* node;
    double possibleDualValue;
    int i;

    assert(modeldata != NULL);
    assert(dualvalues != NULL);
    assert(0 <= day && day < modeldata->nDays);

    possibleDualValue = 0;
    /* if the instance is preprocessed, the customers of this day are exactly the neighbors of the depot */
    if (modeldata->neighbors != NULL)
    {
        node = modeldata->neighbors[modeldata->nC - 1][day];
        while(node != NULL)
        {
            int hardCustomerBonus = 0;
            i = node->id;
            /* Price Collecting for hard customers */
            if (probdata->useOptionals == TRUE && probdata->optionalCustomers[i] == TRUE)
            {
               hardCustomerBonus = modeldata->obj[i] * PRICE_COLLECTING_WEIGHT;
            }
            if (!SCIPisSumPositive(scip, dualvalues[i] + hardCustomerBonus))
            {
                node = node->next;
                continue;
            }
            possibleDualValue += dualvalues[i] + hardCustomerBonus;
            node = node->next;
        }
    /* if the instance is not preprocessed, the customers of this day are searched manually */
    } else {
        for (i = 0; i < modeldata->nC - 1; i++)
        {
           /* Price Collecting for hard customers */
            int hardCustomerBonus = 0;
            if (probdata->useOptionals == TRUE && probdata->optionalCustomers[i] == TRUE)
            {
               hardCustomerBonus = modeldata->obj[i] * PRICE_COLLECTING_WEIGHT;
            }
            if (!SCIPisSumPositive(scip, dualvalues[i] + hardCustomerBonus))
            {
                continue;
            }
            window = modeldata->timeWindows[i];
            while (window != NULL)
            {
                if (window->day == day)
                {
                    possibleDualValue += dualvalues[i] + hardCustomerBonus;
                    break;
                }
                window = window->next;
            }
        }
    }

    assert(possibleDualValue >= 0);
    return possibleDualValue;
}

/**@} */

