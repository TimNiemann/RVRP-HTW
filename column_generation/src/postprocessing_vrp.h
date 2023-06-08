/**@file   postprocessing_vrp_vrp.h
 * @brief  solution structure and functions to create an output file
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __POSTPROCESSING_VRP__
#define __POSTPROCESSING_VRP__

#include "scip/scip.h"
#include "tools_data.h"


SCIP_RETCODE createSolutionWindow(
    SCIP*               scip,
    solutionWindow**    solutionwindow,
    int                 day,
    int                 start,
    int                 end,
    int                 latestArrival,
    modelWindow*        window
);

SCIP_RETCODE freeSolutionWindow(
    SCIP*               scip,
    solutionWindow*     solutionwindow
);

SCIP_RETCODE freeSolutionWindowArray(
    SCIP*               scip,
    solutionWindow**    solutionwindow,
    int                 size
);

SCIP_RETCODE createOutputFile(
    SCIP*               scip,
    model_data*         modeldata,
    char*               outputFile
);


#endif