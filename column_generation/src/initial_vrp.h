/**@file   initial_vrp.h
 * @brief  Methods for generating an initial solution
 * @author Lukas Schürmann, University Bonn
 */

#ifndef __SCIP_INITIAL_VRP__
#define __SCIP_INITIAL_VRP__

#include "scip/scip.h"
#include "tools_data.h"
#include "probdata_vrp.h"


typedef struct _inttuple
{
    int value;
    int index;
} inttuple;

/** improves and adds a set of tours to SCIP */
extern
SCIP_RETCODE addToursToModel(
        SCIP*           scip,
        SCIP_PROBDATA*  probdata,
        int*            dayofnode,
        int**           tours,
        int*            tourlength,
        double*         tourobj,
        char*           algoName,
        SCIP_Bool       callback
);

/** creates and add variables for a given initial solution */
extern
SCIP_RETCODE createSolutionColumns(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata,
        solution_data*        solutiondata
);

/** Adds a set of valid variables to the SPF */
extern
SCIP_RETCODE createGreedyColumns(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata
);

/** Adds one empty route per day */
extern
SCIP_RETCODE createInitialColumns(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata
);

/** dispatching heuristic to find feasible solution */
extern
SCIP_RETCODE initialDispatching(
        SCIP*                 scip,
        SCIP_PROBDATA*        probdata
);

#endif
