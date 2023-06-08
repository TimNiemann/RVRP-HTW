/**@file   label_vrp.h
 * @brief  label data structure and label propagation function
 * @author Tim Niemann, TU Braunschweig
 */


#ifndef __LABEL_VRP_H__
#define __LABEL_VRP_H__

#include "scip/scip.h"
#include "scip/cons_setppc.h"

#include "tools_data.h"

typedef struct _labelVrp
{
    int             node;                   /** current node at which this label is present */
    double          redcost;                /** accumulated reduced costs on this path */
    double          collactableRedcost;     /** reduced costs which could still be collected */
    int*            bitVisitednodes;        /** nC-bit-array of visitednodes */
    int             sizeBitarray;           /** = (modeldata->nC/INT_BIT_SIZE + 1) */
    int*            visitednodes;           /** sequence of all already visited nodes */
    int             nvisitednodes;          /** size of visitednodes */
    int             starttime;              /** time at which the tour started at the depot (currently no further use!) */
    int             day;                    /** day on which this tour is planned */
    int*            arrivaltimes;           /** array of all (delayed) arrival times at this node */ 
    int             narrivaltimes;          /** size of arrivaltimes, equal to robustness parameter + 1 */
    int             ndominated;
    SCIP_Bool       isPropagated;
    SCIP_Bool       isDominated;
    double          lhs;
} labelVrp;

/** Creates an empty label at the depot node */
SCIP_RETCODE labelVrpCreateEmpty(
    SCIP*           scip, 
    labelVrp**      label,
    int             ncustomers,
    int             narrivaltimes,
    double          initialRedcost,
    double          collactableRedcost,
    int             day
);

/** Creates a label at node */
SCIP_RETCODE labelVrpCreate(
    SCIP*           scip, 
    labelVrp**      label,
    int             node,
    int             ncustomers,
    int*            visitednodes,
    int             nvisitednodes,
    int*            arrivaltimes,
    int             narrivaltimes,
    double          redcost,
    double          collactableRedcost,
    int             starttime,
    int             day
);

/** Free label data */
SCIP_RETCODE labelVrpFree(
    SCIP*           scip, 
    labelVrp**      label
);

/** Propagate a label to destination */
SCIP_RETCODE labelVrpPropagate(
    SCIP*           scip,
    model_data*     modeldata,
    labelVrp*       oldLabel,
    labelVrp**      newLabel,
    int             destination,
    double          dualvalue,
    SCIP_Bool       isFarkas
);

/**
 * Checks if the first label dominates the second one
 * @param labelA first label
 * @param labelB second label
 * @return TRUE if labelA dominates labelB, 
 *         FALSE else */
SCIP_Bool labelVrpDominates(
    SCIP*       scip,
    labelVrp*   labelA,
    labelVrp*   labelB
);

/**
 * computes still collectable dualvalues for a given label
 * @param label label with arivaltimes and current day
 * @param dualvalues dualvalues
 * @param upperTimeWindows nC-array with the last possible arrivaltime for each customer, 
 *                         this equals upper boundary of time window when customer has time window on this day
 *                         in all other cases the value is zero
 * @return negative sum over all reduced cost that still could be collected for this label*/
double labelVrpCollactableRedCostTimeDependent(
    SCIP*           scip,
    model_data*     modeldata,
    labelVrp*       label,
    double*         dualvalues,
    int*            upperTimeWindows
);

/** Print label data */
void labelVrpPrint(
    labelVrp*      label
);

#endif