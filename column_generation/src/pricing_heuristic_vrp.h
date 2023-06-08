/**@file   pricing_heurisitc_vrp.h
 * @brief  heuristic pricing as an alternative to the labeling algorithmn
 * @author Tim Niemann, TU Braunschweig
 * @author Lukas Schürmann, University Bonn
 */

#ifndef __PRICING_HEURISTIC_VRP__
#define __PRICING_HEURISTIC_VRP__

#include "scip/scip.h"

/** calls local search pricing heuristic */
SCIP_RETCODE heuristicPricing(
    SCIP*                scip,
    SCIP_Bool            isFarkas        /**< TRUE for farkas-pricing, FALSE for redcost-pricing */    
    );

#endif