/**@file   labeling_algorithm_vrp.h
 * @brief  labeling Algorithm for column generation 
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __LABELING_ALGORITHM_VRP__
#define __LABELING_ALGORITHM_VRP__

#include "scip/scip.h"
#include "labellist_vrp.h"
#include "tools_vrp.h"

/** struct to pass arguments for labeling to worker threads */
typedef struct _arg_struct {
   SCIP*                scip;
   SCIP_Bool            isFarkas;        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
   SCIP_Bool            isHeuristic;
   int                  day;
   SCIP_Bool*           visited;
   SCIP_Bool*           toDepot;
   label_list*          bestLabels;
} arg_struct;

/** labeling algorithm on one thread */
SCIP_RETCODE labelingAlgorithmIterativ(
        SCIP*                scip,
        SCIP_Bool            isFarkas,        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
        SCIP_Bool            isHeuristic,
        int                  nDays,
        tuple*               days,
        SCIP_Bool*           visited,
        SCIP_Bool*           toDepot
   );

/** Same as labeling Algorithm, but uses a different thread for each day */
SCIP_RETCODE labelingAlgorithmParallel(
   SCIP*                scip,
   SCIP_Bool            isFarkas,        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */
   SCIP_Bool            isHeuristic,
   int                  nDays,
   SCIP_Bool*           visited,
   SCIP_Bool*           toDepot
   );

#endif