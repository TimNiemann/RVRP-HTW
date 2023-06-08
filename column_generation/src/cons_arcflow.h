/**@file   cons_arcflow.h
 * @brief  arcflow branching rule for the VRP
 * @author Lukas Schürmann, University Bonn
 */

#ifndef COLUMN_GENERATION_CONSHDLR_H
#define COLUMN_GENERATION_CONSHDLR_H

#include "scip/scip.h"

enum ConsType
{
    PROHIBIT = 0,
    ENFORCE = 1
};
typedef enum ConsType CONSTYPE;

/** creates the arc flow Conshdlr and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrArcFlow(
        SCIP*               scip
        );

/** creates an arc flow constraint */
SCIP_RETCODE SCIPcreateConsArcFlow(
        SCIP*               scip,
        SCIP_CONS**         cons,
        const char*         name,
        int                 tail,
        int                 head,
        CONSTYPE            type,
        SCIP_NODE*          node,
        SCIP_Bool           local
        );

/** returns tail of the arc */
int SCIPgetTailArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** returns head of the arc */
int SCIPgetHeadArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** return constraint type PROHIBIT or ENFORCE */
CONSTYPE SCIPgetTypeArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< samediff constraint */
);

/** returns number of propagated variables for cons */
int SCIPgetNPropagated(
        SCIP*               scip,
        SCIP_CONS*          cons
);

#endif
