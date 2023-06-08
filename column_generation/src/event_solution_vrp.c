/**@file   event_solution_vrp.c
 * @brief  SCIP callback functions to improve a new best solution
 * @author Tim Niemann, TU Braunschweig
 */

#include "event_solution_vrp.h"
#include "scip/scip.h"
#include "tools_data.h"
#include "probdata_vrp.h"
#include "tools_vrp.h"
#include "vardata_vrp.h"
#include "initial_vrp.h"
#include <string.h>

#define EVENTHDLR_NAME         "bestsol"
#define EVENTHDLR_DESC         "event handler for best solutions found"

static
SCIP_RETCODE improveBestSolution(
    SCIP*               scip
    )
{
    SCIP_VAR** vars = NULL;
    SCIP_VARDATA** vardatas = NULL;
    SCIP_SOL* bestsol = NULL;
    double* vals = NULL;
    int nvars, nbinvars;
    int i, j, k;
    model_data* modeldata;
    double* obj = NULL;
    int** tour;
    int* tourLength = NULL;
    int* dayofnode = NULL;
    int numones = 0;
    char name[] = "callbackImprovement";

    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

    if(probdata == NULL)
    {
        SCIPwarningMessage(scip, "Undefined behaviour: probdata = NULL in event_solution.\n");
        return SCIP_OKAY;
    }

    bestsol = SCIPgetBestSol(scip);

    assert(bestsol != NULL);
    modeldata = probdata->modeldata;

    SCIP_CALL( SCIPgetSolVarsData(scip, bestsol, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
    assert(nvars == nbinvars);

    SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &vardatas, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &obj, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tourLength, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tour, modeldata->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &dayofnode, modeldata->nC) );

    SCIP_CALL (SCIPgetSolVals(scip, bestsol, nvars, vars, vals) );

    /* generate tours from solution */
    for(i = 0; i < modeldata->nC; i++)
    {
        dayofnode[i] = -1;
    }
    for (i = 0; i < nvars; i++)
    {
        if (SCIPisEQ(scip, vals[i], 1.0))
        {
            numones++;

            j = SCIPvarGetData(vars[i])->day;
            vardatas[j] = SCIPvarGetData(vars[i]);
            obj[j] = SCIPvarGetObj(vars[i]);
            tourLength[j] = vardatas[j]->tourlength;
            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &tour[j], modeldata->nC) );
            for (k = 0; k < tourLength[j]; k++)
            {
                tour[j][k] = vardatas[j]->customertour[k];

                assert(dayofnode[vardatas[j]->customertour[k]] == -1);

                dayofnode[vardatas[j]->customertour[k]] = j;
                assert(!nodeIsDepot(modeldata, tour[j][k]));
            }
        }
    }
    assert(numones == modeldata->nDays);

    /* try to improve it and add new tours to model */
    SCIP_CALL( addToursToModel( scip, probdata, dayofnode, tour, tourLength, obj, name, TRUE));

    for(i = 0; i < modeldata->nDays; i++)
    {
        SCIPfreeBlockMemoryArray(scip, &tour[i], modeldata->nC);
    }
    SCIPfreeBlockMemoryArray(scip, &dayofnode, modeldata->nC);
    SCIPfreeBlockMemoryArray(scip, &tour, modeldata->nDays);
    SCIPfreeBlockMemoryArray(scip, &tourLength, modeldata->nDays);
    SCIPfreeBlockMemoryArray(scip, &obj, modeldata->nDays);
    SCIPfreeBlockMemoryArray(scip, &vardatas, modeldata->nDays);
    SCIPfreeBufferArray(scip, &vals);

    return SCIP_OKAY;
}


/** copy method for event handler plugins (called when SCIP copies plugins) */
static
SCIP_DECL_EVENTCOPY(eventCopyBestsol) 
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* call inclusion method of event handler */
   SCIP_CALL( SCIPincludeEventHdlrBestsol(scip) );

   return SCIP_OKAY;
}

/** initialization method of event handler (called after problem was transformed) */
static
SCIP_DECL_EVENTINIT(eventInitBestsol)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);

   /* notify SCIP that your event handler wants to react on the event type best solution found */
   SCIP_CALL( SCIPcatchEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** deinitialization method of event handler (called before transformed problem is freed) */
static
SCIP_DECL_EVENTEXIT(eventExitBestsol)
{
   assert(scip != NULL);
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   
   /* notify SCIP that your event handler wants to drop the event type best solution found */
   SCIP_CALL( SCIPdropEvent( scip, SCIP_EVENTTYPE_BESTSOLFOUND, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecBestsol)
{
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(scip != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_BESTSOLFOUND);

   /* only at root */
   if(SCIPnodeGetNumber(SCIPgetCurrentNode(scip)) > 1)
       return SCIP_OKAY;

   SCIPdebugMsg(scip, "exec method of event handler for best solution found\n");
   
   improveBestSolution(scip);
   
   return SCIP_OKAY;
}

/** includes event handler for best solution found */
SCIP_RETCODE SCIPincludeEventHdlrBestsol(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_EVENTHDLRDATA* eventhdlrdata;
   SCIP_EVENTHDLR* eventhdlr;
   eventhdlrdata = NULL;
   
   eventhdlr = NULL;
   /* create event handler for events on watched variables */
   SCIP_CALL( SCIPincludeEventhdlrBasic(scip, &eventhdlr, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecBestsol, eventhdlrdata) );
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPsetEventhdlrCopy(scip, eventhdlr, eventCopyBestsol) );
   SCIP_CALL( SCIPsetEventhdlrInit(scip, eventhdlr, eventInitBestsol) );
   SCIP_CALL( SCIPsetEventhdlrExit(scip, eventhdlr, eventExitBestsol) );
   
   return SCIP_OKAY;
}