/**@file   probdata_vrp.c
 * @brief  problem data for the master problem 
 * @author Tim Niemann, TU Braunschweig
 */

#include <string.h>

#include "probdata_vrp.h"
#include "tools_data.h"
#include "vardata_vrp.h"
#include "pricer_vrp.h"
#include "tools_vrp.h"
#include "initial_vrp.h"

#include "scip/cons_setppc.h"
#include "scip/scip.h"

/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   SCIPdebugMsg(scip, "exec method of event handler for added variable to probdata\n");

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event)) );

   return SCIP_OKAY;
}

/**@} */


/**@name Local methods
 *
 * @{
 */

/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   SCIP_VAR**            vars,               /**< all exist variables */
   SCIP_CONS**           conss,              /**< set partitioning constraints */
   int                   nvars,              /**< number of variables */
   int                   nconss,             /**< number of constraints */
   model_data*           modeldata,          /**< modeldata */
   int                   delayTolerance,     /**< delayTolerance or appointment window length */
   double*               alphas,             /**< weight for balancing the two objective functions */
   SCIP_Bool*            optionalCustomers       /**< nC-Array for which customers are not necessary for a solution */
   )
{
   assert(scip != NULL);
   assert(modeldata != NULL);

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   if( nvars > 0 )
   {
      /* copy variable array */
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->vars, vars, nvars) );
   }
   else
      (*probdata)->vars = NULL;

   /* duplicate arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->conss, conss, nconss) );

   (*probdata)->modeldata = modeldata;
   (*probdata)->nvars = nvars;
   (*probdata)->varssize = nvars;
   (*probdata)->nconss = nconss;
   (*probdata)->delayTolerance = delayTolerance;
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->alphas, 3));
    (*probdata)->alphas[0] = alphas[0];
    (*probdata)->alphas[1] = alphas[1];
    (*probdata)->alphas[2] = alphas[2];
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->shortestEdge, modeldata->nDays));
    for(int i = 0; i < modeldata->nDays; i++)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &((*probdata)->shortestEdge[i]), modeldata->nC));
    }
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->optionalCustomers, optionalCustomers, modeldata->nC -1));
   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;

   assert(scip != NULL);
   assert(probdata != NULL);

   /* release all variables */
   for( i = 0; i < (*probdata)->nvars; ++i )
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->vars[i]) );
   }

   /* release all constraints */
   for( i = 0; i < (*probdata)->nconss; ++i )
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conss[i]) );
   }
   /* free memory of arrays */
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->alphas, 3);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->vars, (*probdata)->varssize);
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->conss, (*probdata)->nconss);
   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigVrp)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransVrp)
{
   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->vars, sourcedata->conss, 
         sourcedata->nvars, sourcedata->nconss, sourcedata->modeldata, sourcedata->delayTolerance, sourcedata->alphas, sourcedata->optionalCustomers) );
   (*targetdata)->optionalCost = sourcedata->optionalCost;
   (*targetdata)->useOptionals = sourcedata->useOptionals;
   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nconss, (*targetdata)->conss, (*targetdata)->conss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, (*targetdata)->nvars, (*targetdata)->vars, (*targetdata)->vars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransVrp)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolVrp)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolVrp)
{  /*lint --e{715}*/
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   model_data*           modeldata,          /**< model data */
   solution_data*        solutiondata,       /**< solution data (optional input) */
   int                   delayTolerance,     /**< delayTolerance or appointment window length */
   double*               alphas,             /**< weight for balancing the two objective functions */
   SCIP_Bool*            optionalCustomers       /**< nC-Array for which customers are not necessary for a solution */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONS** conss;
   char name[SCIP_MAXSTRLEN];
   double optionalCost = 0;
   SCIP_Bool useOptionals = FALSE;
   int i;

   assert(scip != NULL);
   assert(modeldata != NULL);
   assert(modeldata->nDays <= modeldata->nC - 1);
   assert(delayTolerance >= 0);
   SCIPdebugMsg(scip, "number of days for this vehicle routing problem %d", modeldata->nDays);

   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, "Vehicle Routing Problem - Set Partitioning Formulation") );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigVrp) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransVrp) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransVrp) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolVrp) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolVrp) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* allocate memory for the constraints: customer and days */
   SCIP_CALL( SCIPallocBufferArray(scip, &conss, modeldata->nC - 1 + modeldata->nDays) );

   /* create set partitioning constraints for each customer but the depot */
   assert(modeldata->nC > 1);
   for( i = 0; i < modeldata->nC - 1; ++i )
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "customer_%d", modeldata->customerIDs[i]);
      if (SAME_OBJECTIVES)
      {
         modeldata->obj[i] = 1.0;
      }
      if (optionalCustomers[i] == FALSE)
      {
         SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &conss[i], name, 0, NULL) );
      } else {
         SCIP_CALL( SCIPcreateConsBasicSetpack(scip, &conss[i], name, 0, NULL) );
         useOptionals = TRUE;
         optionalCost += modeldata->obj[i];
      }

      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, conss[i], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );   
   }
   assert(i == modeldata->nC -1);
   
   assert(modeldata->nDays > 0);
   /* add constraint, such that exactly one route is choosen for every day */
   for (i = 0; i < modeldata->nDays; i++)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "day_%d", i);
      SCIP_CALL( SCIPcreateConsBasicSetpart(scip, &conss[modeldata->nC - 1 + i], name, 0, NULL) );
      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, conss[modeldata->nC - 1 + i], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, conss[modeldata->nC - 1 + i]) );   
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, NULL, conss, 0, modeldata->nC - 1 + modeldata->nDays, modeldata, delayTolerance, alphas, optionalCustomers) );
   assert(probdata != NULL);

   /* set recovery data */
   probdata->optionalCost = optionalCost;
   probdata->useOptionals = useOptionals;
   if (useOptionals == FALSE)
   {
      probdata->optionalCost = SCIP_DEFAULT_INFINITY;
   }

   /* set user problem data and pricer */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );
   SCIP_CALL( SCIPpricerVrpActivate(scip, conss, probdata->nconss, modeldata) );

   /* create initial columns, in the best case create an initial solution */
   SCIP_CALL(createInitialColumns(scip, probdata));
   if(solutiondata != NULL)
   {
       printf("Add solution input as initial columns\n");
       SCIP_CALL( createSolutionColumns(scip, probdata, solutiondata) );
   }else {
       printf("Find initial tours:\n");

       SCIP_CALL(createGreedyColumns(scip, probdata));

       SCIP_CALL(initialDispatching(scip, probdata));
   }

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &conss);
   
   return SCIP_OKAY;
}


/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   /* check if enough memory is left */
   if( probdata->varssize == probdata->nvars )
   {
      int newsize;
      newsize = MAX(100, probdata->varssize * 2);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->vars, probdata->varssize, newsize) );
      probdata->varssize = newsize;
   }

   /* caputure variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   probdata->vars[probdata->nvars] = var;
   probdata->nvars++;

   SCIPdebugMsg(scip, "added variable to probdata; nvars = %d\n", probdata->nvars);

   return SCIP_OKAY;
}

/**
 * Check by name, if a potential new variable was already added to the master problem
 * @param name variable name to check
 * @param probdata problem data which contains all previously generated variables
 * @return TRUE, if this variable is already contained in the master problem;
 *         FALSE, if this varibale name is not found */
extern
SCIP_Bool SCIPprobdataContainsVar(
   SCIP_PROBDATA*        probdata,           /**< problem data */
   const char*           name                /**< variables name to check */
   )
{
   int i;
   assert(probdata != NULL);
   assert(name != NULL);

   for (i = 0; i < probdata->nvars; i++)
   {
      const char* tempname = SCIPvarGetName(probdata->vars[i]);
      if (strcmp(name, tempname) == 0)
      {
         SCIPdebugMessage("Tried to add variable that is already contained in master problem: %s\n", name);
         return TRUE;
      }
   }
   return FALSE;
}

/**@} */
