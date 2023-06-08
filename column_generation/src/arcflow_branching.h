/**@file   branching.h
 * @brief  branching rule for the VRP
 * @author Lukas Schürmann, University Bonn
 */

#ifndef __ARCFLOW_BRANCHING__
#define __ARCFLOW_BRANCHING__

#include "scip/scip.h"
#include "tools_data.h"

/** creates the arc flow branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleArcFlow(
        SCIP*                 scip                /**< SCIP data structure */
);

#endif
