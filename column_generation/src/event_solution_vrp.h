/**@file   event_solution_vrp.h
 * @brief  SCIP callback functions to improve a new solution
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __EVENT_SOLUTION_VRP__
#define __EVENT_SOLUTION_VRP__

#include "scip/scip.h"

/** includes event handler for new solution found */
SCIP_RETCODE SCIPincludeEventHdlrBestsol(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif