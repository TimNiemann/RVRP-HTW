/**@file   vardata_vrp.c
 * @brief  varibale data for the variables of the master problem, each represents one tour
 * @author Tim Niemann, TU Braunschweig
 */

#include <string.h>
#include "scip/scip.h"
#include "scip/cons_setppc.h"
#include "probdata_vrp.h"
#include "vardata_vrp.h"
#include "tools_vrp.h"
#include "postprocessing_vrp.h"


/**@name Local methods
 *
 * @{
 */

/** create a vardata */
static
SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  customertour,       /**< sequence of visited customers on this tour (without depot at start/end) */ 
   int                   tourlength,         /**< number of visited customers on this tour (without depot) */
   int                   tourduration,       /**< worst case cummulated time from start at depot until arrival back at the depot in seconds */
   solutionWindow**      solutionwindows,
   int                   day                 /**< number of days after beginning of planning horizon for when this tour is planned */
   )
{
   SCIP_CALL( SCIPallocMemory(scip, vardata) );

   assert(customertour != NULL);
   assert(tourlength >= 0);
   assert(tourduration >= 0);

   if (tourlength > 0)
   {
      int i;
      SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*vardata)->customertour, customertour, tourlength) );
      SCIP_CALL( SCIPallocMemoryArray(scip, &(*vardata)->solutionwindows, tourlength) );
      for (i = 0; i < tourlength; i++)
      {
         int start = 0;
         int end = 0;
         int latestArrival = 0;
         modelWindow* window = NULL;
         if (solutionwindows != NULL)
         {
            start = solutionwindows[i]->start_t;
            end = solutionwindows[i]->end_t;
            latestArrival = solutionwindows[i]->latest_arrival;
            window = solutionwindows[i]->modelwindow;
         }
         SCIP_CALL( createSolutionWindow(scip, &(*vardata)->solutionwindows[i], day, start, end, latestArrival, window));
      }
   } else {
      (*vardata)->customertour = NULL;
      (*vardata)->solutionwindows = NULL;
   }

   (*vardata)->tourlength = tourlength;
   (*vardata)->tourduration = tourduration;
   (*vardata)->day = day;

   return SCIP_OKAY;
}

/** frees user data of variable */
static
SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   if ((*vardata)->customertour != NULL)
   {
      SCIPfreeMemoryArray(scip, &(*vardata)->customertour);
   }
   if ((*vardata)->solutionwindows != NULL)
   {
      SCIP_CALL( freeSolutionWindowArray(scip, (*vardata)->solutionwindows, (*vardata)->tourlength));
   }
   SCIPfreeMemory(scip, vardata);
   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees user data of transformed variable (called when the transformed variable is freed) */
static
SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** create variable data */
SCIP_RETCODE SCIPvardataCreateVrp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  customertour,       /**< sequence of visited customers on this tour (without depot at start/end) */ 
   int                   tourlength,         /**< number of visited customers on this tour (without depot) */
   int                   tourduration,       /**< worst case cummulated time from start at depot until arrival back at the depot in seconds */
   solutionWindow**      solutionwindows,
   int                   day                 /**< number of days after beginning of planning horizon for when this tour is planned */
   )
{
   SCIP_CALL( vardataCreate(scip, vardata, customertour, tourlength, tourduration, solutionwindows, day) );

   return SCIP_OKAY;
}

/** creates variable */
SCIP_RETCODE SCIPcreateVarVrp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   /* create a basic variable object */
   SCIP_CALL( SCIPcreateVarBasic(scip, var, name, 0.0, 1.0, obj, SCIP_VARTYPE_BINARY) );
   assert(*var != NULL);

   /* set callback functions */
   SCIPvarSetData(*var, vardata);
   SCIPvarSetDeltransData(*var, vardataDelTrans);

   /* set initial and removable flag */
   SCIP_CALL( SCIPvarSetInitial(*var, initial) );
   SCIP_CALL( SCIPvarSetRemovable(*var, removable) );

   SCIPvarMarkDeletable(*var);

   return SCIP_OKAY;
}

SCIP_RETCODE SCIPcreateColumn(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata,
        const char*           name,
        SCIP_Bool             initial,
        double                obj,
        int*                  tour,
        int                   tourlength,
        int                   tourduration,
        solutionWindow**      solutionwindows,
        int                   day
){
    SCIP_CONS** conss;
    SCIP_VARDATA* vardata;
    SCIP_VAR* var;
    int i;

    assert(probdata != NULL);
    assert(probdata->modeldata != NULL);
    assert(probdata->conss != NULL);

    assert(tour != NULL);
    assert(solutionwindows != NULL || tourlength == 0);
    assert(tourlength >= 0);
    assert(day >= 0);
    assert(tourduration <= WORKTIME_LIMIT);

    conss = probdata->conss;

    if (PRINT_GENERATED_TOURS) {
        printf("%s: %f\n", name, obj);
    }
    //Use the tour for both arrays although not the whole array is copied for the tour: just the first entry
    SCIP_CALL( SCIPvardataCreateVrp(scip, &vardata, tour, tourlength, tourduration, solutionwindows, day ) );

    /* create variable for the vrp which contains only this customer */
    SCIP_CALL( SCIPcreateVarVrp(scip, &var, name, obj, initial, TRUE, vardata) );

    /* add variable to the problem */
    if (initial)
    {
        SCIP_CALL( SCIPaddVar(scip, var) );
    } else {
        SCIP_CALL( SCIPaddPricedVar(scip, var, 1.0) );
    }

    /* change the upper bound of the binary variable to lazy since the upper bound is already enforced
     * due to the objective function in the set partitioning constraint;
     * The reason for doing is that, is to avoid the bound of x <= 1 in the LP relaxation since this bound
     * constraint would produce a dual variable which might have a positive reduced cost
     */
    SCIP_CALL( SCIPchgVarUbLazy(scip, var, 1.0) );

    /* store variable in the problme data */
    if (initial)
    {
        SCIP_CALL( SCIPprobdataAddVar(scip, probdata, var) );
    }

    for(i = 0; i < tourlength; i++)
    {
        assert(0 <= tour[i] && tour[i] < probdata->modeldata->nC - 1);
        assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(conss[tour[i]]) ), "setppc", 6) );
        SCIP_CALL(SCIPaddCoefSetppc(scip, conss[tour[i]], var));
    }
    /* add variable to corresponding set partitioning constraint */
    assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(conss[probdata->modeldata->nC - 1 + day]) ), "setppc", 6) );   
    SCIP_CALL( SCIPaddCoefSetppc(scip, conss[probdata->modeldata->nC - 1 + day], var) );

    /* add the variable data to the variable */
    SCIPvarSetData(var, vardata);



    /* release variable */
    SCIP_CALL( SCIPreleaseVar(scip, &var) );

    return SCIP_OKAY;
}