/**@file   probdata_vrp.h
 * @brief  problem data for the master problem 
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __SCIP_PROBDATA_VRP__
#define __SCIP_PROBDATA_VRP__

#include "scip/scip.h"
#include "tools_data.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the vrp, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   SCIP_VAR**            vars;               /**< all exiting variables in the problem */
   SCIP_CONS**           conss;              /**< set partitioning constraints for each customer exactly one and one for the number of tours*/
   int                   nvars;              /**< number of generated variables */
   int                   varssize;           /**< size of the variable array, always greater or equals to nvars */
   int                   nconss;             /**< number of constraints */
   model_data*           modeldata;          /**< model data */
   int                   delayTolerance;     /**< delayTolerance or appointment window length */
   double                weightObjective;    /**< weight for balancing the two objective functions */
   double*               alphas;
   int**                 shortestEdge;
   SCIP_Bool             useOptionals;       /**< if set, the optional customers won't be necessary */
   SCIP_Bool*            optionalCustomers;  /**< nC-Array for which customers are not necessary for a solution */
   double                optionalCost;       /**< sum of all optional customer objective values */
};

/** sets up the problem data */
extern
SCIP_RETCODE SCIPprobdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   model_data*           modeldata,          /**< model data */
   solution_data*        solutiondata,       /**< solution data (optional input) */
   int                   delayTolerance,     /**< delayTolerance or appointment window length */
   double*               alphas,             /**< weight for balancing the two objective functions */
   SCIP_Bool*            optionalCustomers       /**< nC-Array for which customers are not necessary for a solution */
   );

/** adds given variable to the problem data */
extern
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   );

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
   );
#endif
