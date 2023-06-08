/**@file   primal_heuristic_vrp.c
 * @brief  Functions for generating an inital solution
 * @author Lukas Schürmann, University Bonn
 */

#ifndef __PRIMAL_HEURISTIC_VRP__
#define __PRIMAL_HEURISTIC_VRP__

#include "scip/scip.h"

/** heuristic to calculate a good VRP solution based current LP-solution */
SCIP_RETCODE dispatchingHeuristic(
        SCIP*               scip,
        SCIP_PROBDATA*      probdata,
        SCIP_PRICERDATA*    pricerdata
);

#endif