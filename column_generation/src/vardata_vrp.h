/**@file   vardata_vrp.h
 * @brief  varibale data for the variables of the master problem, each represents one tour
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __SCIP_VARDATA_VRP__
#define __SCIP_VARDATA_VRP__

#include "scip/scip.h"
#include "postprocessing_vrp.h"

/** Variable data which is attached to all variables.
 *
 *  This variable data is used to store in which constraints this variable appears. Therefore, the variable data
 *  contains the ids of constraints in which the variable is part of. Hence, that data give us a column view.
 */
struct SCIP_VarData
{
   int                   day;                /**< number of days after beginning of planning horizon for when this tour is planned, -1 if not specified */
   int*                  customertour;       /**< sequence of visited customers on this tour (without depot at start/end) */ 
   int                   tourlength;         /**< number of visited customers on this tour (without depot) */
   int                   tourduration;       /**< worst case cummulated time from start at depot until arrival back at the depot in seconds */
   solutionWindow**      solutionwindows;
};

/** create variable data */
extern
SCIP_RETCODE SCIPvardataCreateVrp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  customertour,       /**< sequence of visited customers on this tour (without depot at start/end) */ 
   int                   tourlength,         /**< number of visited customers on this tour (without depot) */
   int                   tourduration,       /**< worst case cummulated time from start at depot until arrival back at the depot in seconds */
   solutionWindow**      solutionwindows,
   int                   day                 /**< number of days after beginning of planning horizon for when this tour is planned, -1 if not specified */
   );

/** creates variable */
extern
SCIP_RETCODE SCIPcreateVarVrp(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** Adds a new variable to the scip model for a given tour */
extern
SCIP_RETCODE SCIPcreateColumn(
   SCIP*                 scip,              /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,          /**< Probdata structer */
   const char*           name,              /**< name of variable */
   SCIP_Bool             initial,           /**< should var's column be present in the initial root LP?*/
   double                obj,               /**< objective function value */
   int*                  tour,              /**< array of visited customers */
   int                   tourlength,        /**< length of the tour */
   int                   tourduration,      /**< tour duration in seconds */
   solutionWindow**      solutionwindows,
   int                   day                /**< day of the tour */
);

#endif
