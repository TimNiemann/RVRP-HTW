/**@file   cons_arcflow.c
 * @brief  arcflow branching rule for the VRP
 * @author Lukas Schürmann, University Bonn
 */
#include <assert.h>

#include "cons_arcflow.h"
#include "probdata_vrp.h"
#include "vardata_vrp.h"

#define CONSHDLR_NAME          "arcflow"
#define CONSHDLR_DESC          "stores the local branching decisions"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/**@} */

/*
 * Data structures
 */

/** Constraint data for  \ref cons_arcflow.c "ArcFlow" constraints */
struct SCIP_ConsData
{
    int                   tail;               /**< tail of the arc */
    int                   head;               /**< head of the arc */
    CONSTYPE              type;               /**< stores whether arc gets enforced or prohibited */
    int                   npropagatedvars;    /**< number of variables that existed, the last time, the related node was
                                              *   propagated, used to determine whether the constraint should be
                                              *   repropagated*/
    int                   npropagations;      /**< stores the number propagations runs of this constraint */
    unsigned int          propagated:1;       /**< is constraint already propagated? */
    SCIP_NODE*            node;               /**< the node in the B&B-tree at which the cons is sticking */
};

/**@name Local methods
 *
 * @{
 */

/** create constraint data */
static
SCIP_RETCODE consdataCreate(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
        int                   tail,               /**< tail of the arc */
        int                   head,               /**< head of the arc */
        CONSTYPE              type,               /**< stores whether arc gets enforced or prohibited */
        SCIP_NODE*            node                /**< the node in the B&B-tree at which the cons is sticking */
)
{
    assert( scip != NULL );
    assert( consdata != NULL );
    assert( tail >= 0 );
    assert( head >= 0 );
    assert( type == ENFORCE || type == PROHIBIT );


    SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

    (*consdata)->tail = tail;
    (*consdata)->head = head;
    (*consdata)->type = type;
    (*consdata)->npropagatedvars = 0;
    (*consdata)->npropagations = 0;
    (*consdata)->propagated = FALSE;
    (*consdata)->node = node;

    return SCIP_OKAY;
}

/** display constraints */
static
void consdataPrint(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONSDATA*        consdata,           /**< constraint data */
        FILE*                 file                /**< file stream */
)
{
    SCIP_PROBDATA* probdata;

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    SCIPinfoMessage(scip, file, "%s(%d,%d) at node %d\n",
                    consdata->type == PROHIBIT ? "prohibit" : "enforce",
                    consdata->tail, consdata->head, SCIPnodeGetNumber(consdata->node) );
}

static
int arcInVar(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_VAR*             var,                /**< variable to check */
        int                   tail,               /**< tail of arc to check */
        int                   head                /**< head of arc to check */
)
{
    int j;
    SCIP_VARDATA* vardata;
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);

    vardata = SCIPvarGetData(var);

    if(vardata->tourlength <= 0)
    {
        return 0;
    }

    /* special case head/tail = depot */
    if(head == probdata->modeldata->nC - 1)
    {
        if(vardata->customertour[vardata->tourlength - 1] == tail)
        {
            return 2;
        }
    }
    if(tail == probdata->modeldata->nC - 1)
    {
        if(vardata->customertour[0] == head)
        {
            return 2;
        }
    }

    /** Possible cases for tail = u, head = v and w != u, w != v:
     *  1. u, v both not in tour
     *  2. tour = ... -> u -> w -> ...
     *  3. tour = ... -> w -> v -> ...
     *  4. tour = ... -> u -> v -> ...
     **/
    if(vardata->customertour[vardata->tourlength - 1] == tail) // case 2
    {
        return 1;
    }
    for(j = 0; j < vardata->tourlength - 1; j++)
    {
        if(vardata->customertour[j] == head) // case 2
        {
            return 1;
        }
        if(vardata->customertour[j] == tail)
        {
            if(vardata->customertour[j+1] != head) // case 3
            {
                return 1;
            }
            // case 4
            return 2;
        }
    }
    if(vardata->customertour[vardata->tourlength - 1] == head) // case 3
    {
        return 1;
    }
    return 0; // case 1
}

/** fixes a variable to zero if the corresponding packings are not valid for this constraint/node (due to branching) */
static
SCIP_RETCODE checkVariable(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONSDATA*        consdata,           /**< constraint data */
        SCIP_VAR*             var,                /**< variables to check  */
        int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
        SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
)
{
    int arcStatus;
    CONSTYPE type;

    SCIP_Bool fixed;
    SCIP_Bool infeasible;

    assert(scip != NULL);
    assert(consdata != NULL);
    assert(var != NULL);
    assert(nfixedvars != NULL);
    assert(cutoff != NULL);

    /* if variables is locally fixed to zero continue */
    if( SCIPvarGetUbLocal(var) < 0.5 )
        return SCIP_OKAY;

    /* check if the tour which corresponds to the variable is feasible for this constraint */
    arcStatus = arcInVar(scip, var, consdata->tail, consdata->head);

    /** arcStatus = 0 (tail and head both not visited by tour) -> variable allowed
     *  arcStatus = 1 (at least on of them is visited but not the arc itself) -> variable prohibited if arc gets enforced
     *  arcStatus = 2 (arc is used) -> variable prohibited if arc gets prohibited  */

    if(arcStatus == 0) return SCIP_OKAY;
    type = consdata->type;
    if( (type == PROHIBIT && arcStatus == 2) || (type == ENFORCE && arcStatus == 1) )
    {
        SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );

        if( infeasible )
        {
            assert( SCIPvarGetLbLocal(var) > 0.5 );
            SCIPdebugMsg(scip, "-> cutoff\n");
            (*cutoff) = TRUE;
        }
        else
        {
            assert(fixed);
            (*nfixedvars)++;
        }
    }

    return SCIP_OKAY;
}

/** fixes variables to zero if the corresponding tours are not valid for this constraint/node (due to branching) */
static
SCIP_RETCODE consdataFixVariables(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONSDATA*        consdata,           /**< constraint data */
        SCIP_VAR**            vars,               /**< generated variables */
        int                   nvars,              /**< number of generated variables */
        SCIP_RESULT*          result              /**< pointer to store the result of the fixing */
)
{
    int nfixedvars;
    int v;
    SCIP_Bool cutoff;

    nfixedvars = 0;
    cutoff = FALSE;

    SCIPdebugMsg(scip, "check variables %d to %d\n", consdata->npropagatedvars, nvars);

    for( v = consdata->npropagatedvars; v < nvars && !cutoff; ++v )
    {
        SCIP_CALL( checkVariable(scip, consdata, vars[v], &nfixedvars, &cutoff) );
    }

    SCIPdebugMsg(scip, "fixed %d variables locally\n", nfixedvars);

    if( cutoff )
        *result = SCIP_CUTOFF;
    else if( nfixedvars > 0 )
        *result = SCIP_REDUCEDDOM;

    return SCIP_OKAY;
}

/** check if all variables are valid for the given consdata */
#ifndef NDEBUG
static
SCIP_Bool consdataCheck(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_PROBDATA*        probdata,           /**< problem data */
        SCIP_CONSDATA*        consdata,           /**< constraint data */
        SCIP_Bool             beforeprop          /**< is this check performed before propagation? */
)
{
    SCIP_VAR** vars;
    int nvars;

    SCIP_VAR* var;

    int arcStatus;
    CONSTYPE type;

    int v;

    vars = probdata->vars;
    nvars = (beforeprop ? consdata->npropagatedvars : probdata->nvars);
    assert(nvars <= probdata->nvars);

    /* search for variables that stand in conflict with cons */
    for( v = 0; v < nvars; ++v )
    {
        var = vars[v];

        /* if variable is locally fixed to zero continue */
        if( SCIPvarGetUbLocal(var) < 0.5 )
            continue;

        arcStatus = arcInVar(scip, var, consdata->tail, consdata->head);

        if(arcStatus == 0) return SCIP_OKAY;

        type = consdata->type;
        if( (type == PROHIBIT && arcStatus == 2) || (type == ENFORCE && arcStatus == 1) )
        {
            SCIPdebug( consdataPrint(scip, consdata, NULL) );
            SCIPdebug( SCIPprintVar(scip, var, NULL) );
            return FALSE;
        }
    }

    return TRUE;
}
#endif

/** frees arcflow constraint data */
static
SCIP_RETCODE consdataFree(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
)
{
    assert(consdata != NULL);
    assert(*consdata != NULL);

    SCIPfreeBlockMemory(scip, consdata);

    return SCIP_OKAY;
}
/**@} */

/**@name Callback methods
 *
 * @{
 */

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteArcflow)
{  /*lint --e{715}*/
            assert(conshdlr != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(consdata != NULL);
    assert(*consdata != NULL);

    /* free arcflow constraint */
    SCIP_CALL( consdataFree(scip, consdata) );

    return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransArcflow)
{
    SCIP_CONSDATA* sourcedata;
    SCIP_CONSDATA* targetdata;

    assert(conshdlr != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
    assert(sourcecons != NULL);
    assert(targetcons != NULL);

    sourcedata = SCIPconsGetData(sourcecons);
    assert(sourcedata != NULL);

    /* create constraint data for target constraint */
    SCIP_CALL( consdataCreate(scip, &targetdata,
    sourcedata->tail, sourcedata->head, sourcedata->type, sourcedata->node) );

    /* create target constraint */
    SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
    SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
    SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
    SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
    SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

    return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
#define consEnfolpArcflow NULL

/** constraint enforcing method of constraint handler for pseudo solutions */
#define consEnfopsArcflow NULL

/** feasibility check method of constraint handler for integral solutions */
#define consCheckArcflow NULL

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropArcflow)
{
    SCIP_PROBDATA* probdata;
    SCIP_CONSDATA* consdata;

    SCIP_VAR** vars;
    int nvars;
    int c;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(result != NULL);

    SCIPdebugMsg(scip, "propagation constraints of constraint handler <"CONSHDLR_NAME">\n");

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    vars = probdata->vars;
    nvars = probdata->nvars;

    *result = SCIP_DIDNOTFIND;
    for( c = 0; c < nconss; ++c )
    {
        consdata = SCIPconsGetData(conss[c]);

        /* check if all previously generated variables are valid for this constraint */
        assert( consdataCheck(scip, probdata, consdata, TRUE) );

#ifndef NDEBUG
        {
            /* check if there are no equal consdatas */
            SCIP_CONSDATA* consdata2;
            int i;

            for( i = c+1; i < nconss; ++i )
            {
                consdata2 = SCIPconsGetData(conss[i]);

                assert( !(consdata->tail != probdata->modeldata->nC - 1 && consdata->tail == consdata2->tail && consdata->type == ENFORCE && consdata2->type == ENFORCE));
                assert( !(consdata->head != probdata->modeldata->nC - 1 && consdata->head == consdata2->head && consdata->type == ENFORCE && consdata2->type == ENFORCE));
                assert( !(consdata->tail == consdata2->tail && consdata->head == consdata2->head) );
            }
        }
#endif
        if( !consdata->propagated )
        {
            SCIPdebugMsg(scip, "propagate constraint <%s> ", SCIPconsGetName(conss[c]));
            SCIPdebug( consdataPrint(scip, consdata, NULL) );
            SCIP_CALL( consdataFixVariables(scip, consdata, vars, nvars, result) );
            consdata->npropagations++;

            if( *result != SCIP_CUTOFF )
            {
                consdata->propagated = TRUE;
                consdata->npropagatedvars = nvars;
            }
            else
                break;
        }

        /* check if constraint is completely propagated */
        assert( consdataCheck(scip, probdata, consdata, FALSE) );
    }

    return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
#define consLockArcflow NULL

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveArcflow)
{
    SCIP_CONSDATA* consdata;
    SCIP_PROBDATA* probdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->npropagatedvars <= probdata->nvars);

    SCIPdebugMsg(scip, "activate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
    SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
    SCIPdebug( consdataPrint(scip, consdata, NULL) );

    if( consdata->npropagatedvars != probdata->nvars )
    {
        SCIPdebugMsg(scip, "-> mark constraint to be repropagated\n");
        consdata->propagated = FALSE;
        SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
    }

    return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveArcflow)
{
    SCIP_CONSDATA* consdata;
    SCIP_PROBDATA* probdata;

    assert(scip != NULL);
    assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);
    assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

    probdata = SCIPgetProbData(scip);
    assert(probdata != NULL);

    SCIPdebugMsg(scip, "deactivate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
    SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
    SCIPdebug( consdataPrint(scip, consdata, NULL) );

    /* set the number of propagated variables to current number of variables is SCIP */
    consdata->npropagatedvars = probdata->nvars;

    return SCIP_OKAY;
}

/** constraint display method of constraint handler */
static
SCIP_DECL_CONSPRINT(consPrintArcflow)
{
    SCIP_CONSDATA*  consdata;

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    consdataPrint(scip, consdata, file);

    return SCIP_OKAY;
}

/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the arc flow Conshdlr and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrArcFlow(
        SCIP*              scip             /**< SCIP data structure */
)
{
    SCIP_CONSHDLRDATA* conshdlrdata = NULL;
    SCIP_CONSHDLR* conshdlr = NULL;

    /* include constraint handler */
    SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
                                        CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
                                        consEnfolpArcflow, consEnfopsArcflow, consCheckArcflow, consLockArcflow,
                                        conshdlrdata) );
    assert(conshdlr != NULL);

    SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteArcflow) );
    SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransArcflow) );
    SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropArcflow, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
                                   CONSHDLR_PROP_TIMING) );
    SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveArcflow) );
    SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveArcflow) );
    SCIP_CALL( SCIPsetConshdlrPrint(scip, conshdlr, consPrintArcflow) );

    return SCIP_OKAY;
}

/** creates an arc flow constraint */
SCIP_RETCODE SCIPcreateConsArcFlow(
        SCIP*               scip,                /**< SCIP data structure */
        SCIP_CONS**         cons,                /**< pointer to hold the created constraint */
        const char*         name,                /**< name of the constraint */
        int                 tail,                /**< tail of the arc */
        int                 head,                /**< head of the arc */
        CONSTYPE            type,                /**< stores whether arc gets enforced or prohibited */
        SCIP_NODE*          node,                /**< the node in the B&B-tree at which the cons is sticking */
        SCIP_Bool           local                /**< is constraint only valid locally? */
){
    SCIP_CONSHDLR* conshdlr;
    SCIP_CONSDATA* consdata;

    /* find the arcflow constraint handler */
    conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
    if( conshdlr == NULL )
    {
        SCIPerrorMessage("arcflow constraint handler not found\n");
        return SCIP_PLUGINNOTFOUND;
    }

    /* create the constraint specific data */
    SCIP_CALL( consdataCreate(scip, &consdata, tail, head, type, node) );

    /* create constraint */
    SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
                              local, FALSE, FALSE, FALSE, TRUE) );

    SCIPdebugMsg(scip, "created constraint: ");
    SCIPdebug( consdataPrint(scip, consdata, NULL) );

    return SCIP_OKAY;
}

/** returns tail of the arc */
int SCIPgetTailArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< arcflow constraint */
)
{
    SCIP_CONSDATA* consdata;

    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    return consdata->tail;
}

/** returns head of the arc */
int SCIPgetHeadArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< arcflow constraint */
)
{
    SCIP_CONSDATA* consdata;

    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    return consdata->head;
}

/** return constraint type PROHIBIT or ENFORCE */
CONSTYPE SCIPgetTypeArcflow(
        SCIP*                 scip,               /**< SCIP data structure */
        SCIP_CONS*            cons                /**< arcflow constraint */
)
{
    SCIP_CONSDATA* consdata;

    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    return consdata->type;
}

/** returns number of propagated variables for cons */
int SCIPgetNPropagated(
        SCIP*               scip,
        SCIP_CONS*          cons
        ){
    SCIP_CONSDATA* consdata;

    assert(cons != NULL);

    consdata = SCIPconsGetData(cons);
    assert(consdata != NULL);

    return consdata->npropagatedvars;
}

/**@} */
