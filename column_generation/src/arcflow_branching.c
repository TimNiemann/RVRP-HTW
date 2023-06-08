/**@file   branching.c
 * @brief  branching rule for the VRP
 * @author Lukas Schürmann, University Bonn
 */

#include <stdio.h>
#include <assert.h>
#include <math.h>

#include "scip/scip.h"

#include "arcflow_branching.h"
#include "tools_vrp.h"
#include "tools_data.h"
#include "probdata_vrp.h"
#include "vardata_vrp.h"
#include "cons_arcflow.h"

#define BRANCHRULE_NAME            "ArcFlowBranching"
#define BRANCHRULE_DESC            "Branching rule for the arc flow variables"
#define BRANCHRULE_PRIORITY        50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/** check for feasibility of the current node */
static
SCIP_RETCODE feasibilityCheck(
        SCIP*   scip
)
{
    SCIP_VAR** vars;
    SCIP_VAR* var;
    SCIP_VARDATA* vardata;
    SCIP_CONS** conss;
    SCIP_CONS*  cons;
    SCIP_CONSHDLR* conshdlr;
    SCIP_PROBDATA* probdata;
    int nvars, ncons;
    int i, j, k;
    int tail, head;
    CONSTYPE type;
    int arcstatus;

    probdata = SCIPgetProbData(scip);

    conshdlr = SCIPfindConshdlr(scip, "arcflow");

    conss = SCIPconshdlrGetConss(conshdlr);
    ncons = SCIPconshdlrGetNConss(conshdlr);

    SCIP_CALL( SCIPgetVarsData(scip, &vars, &nvars, NULL, NULL, NULL, NULL));

    for(i = 0; i < ncons; i++)
    {
        cons = conss[i];

        if(!SCIPconsIsActive(cons)) continue;

        tail = SCIPgetTailArcflow(scip, cons);
        head = SCIPgetHeadArcflow(scip, cons);
        type = SCIPgetTypeArcflow(scip, cons);

        for(k = 0; k < nvars; k++)
        {
            arcstatus = 0;
            var = vars[k];
            if(SCIPvarGetUbLocal(var) < 0.5) continue;

            vardata = SCIPvarGetData(var);

            if(vardata->tourlength <= 0)
            {
                continue;
            }

            /* special case head/tail = depot */
            if(head == probdata->modeldata->nC - 1)
            {
                if(vardata->customertour[vardata->tourlength - 1] == tail)
                {
                    arcstatus = 2;
                }
            }
            if(tail == probdata->modeldata->nC - 1)
            {
                if(vardata->customertour[0] == head)
                {
                    arcstatus = 2;
                }
            }
            if(arcstatus == 0)
            {
                if(vardata->customertour[vardata->tourlength - 1] == tail) // case 2
                {
                    arcstatus = 1;
                }
                for(j = 0; j < vardata->tourlength - 1; j++)
                {
                    if(vardata->customertour[j] == head) // case 2
                    {
                        arcstatus = 1;
                        break;
                    }
                    if(vardata->customertour[j] == tail)
                    {
                        if(vardata->customertour[j+1] != head) // case 3
                        {
                            arcstatus = 1;
                            break;
                        }
                        // case 4
                        arcstatus = 2;
                        break;
                    }
                }
                if(arcstatus == 0 && vardata->customertour[vardata->tourlength - 1] == head) // case 3
                {
                    arcstatus = 1;
                }
            }

            if(arcstatus == 0) continue;

            if( (type == PROHIBIT && arcstatus == 2) || (type == ENFORCE && arcstatus == 1) )
            {
                printf("Arc (%d, %d), type %s, status %d, VAR %s\n", tail, head, type == PROHIBIT ? "prohibit" : "enforce", arcstatus, SCIPvarGetName(var));
                assert(FALSE);
            }
        }
    }

    return SCIP_OKAY;
}

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpArcFlow)
{
    SCIP_PROBDATA* probdata;
    model_data* modeldata;
    SCIP_Real** arcweights;
    SCIP_VAR** modelvars;
    SCIP_VARDATA* vardata;
    SCIP_Real lpval;
    int nvars;
    int i, j;
    int start, end;
    SCIP_Real bestval = 0.5;

    SCIP_NODE* childprohibit;
    SCIP_NODE* childenforce;
    SCIP_CONS* consprohibit;
    SCIP_CONS* consenforce;

//    SCIP_CALL(feasibilityCheck(scip));

    probdata = SCIPgetProbData(scip);
    modeldata = probdata->modeldata;

    /* allocate arc flow data */
    SCIP_CALL( SCIPallocBlockMemoryArray( scip, &arcweights, probdata->modeldata->nC) );
    for(i = 0; i < probdata->modeldata->nC; i++)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray( scip, &arcweights[i], probdata->modeldata->nC) );
        for(j = 0; j < modeldata->nC; j++)
        {
            arcweights[i][j] = 0.0;
        }
    }

    SCIP_CALL( SCIPgetVarsData(scip, &modelvars, &nvars, NULL, NULL, NULL, NULL));

    /* calculate arc flow values */
    for(i = 0; i < nvars; i++)
    {
        lpval = SCIPvarGetLPSol(modelvars[i]);
        if(SCIPisSumPositive(scip, lpval)
        && SCIPisSumPositive(scip, 1 - lpval))
        {
            vardata = SCIPvarGetData(modelvars[i]);
            if(vardata->tourlength == 0) continue;

            arcweights[modeldata->nC - 1][vardata->customertour[0]] += lpval;
            arcweights[vardata->customertour[vardata->tourlength - 1]][modeldata->nC - 1] += lpval;
            for(j = 1; j < vardata->tourlength; j++)
            {
                arcweights[vardata->customertour[j - 1]][vardata->customertour[j]] += lpval;
            }
        }
    }
    /* search for most fractional arc flow variable */
    for(i = 0; i < modeldata->nC; i++)
    {
        for(j = 0; j < modeldata->nC; j++)
        {
            if(SCIPisSumPositive(scip, arcweights[i][j]) && SCIPisSumPositive(scip, 1 - arcweights[i][j]))
            {
                if(SCIPisSumPositive(scip, bestval - fabs(0.5 - arcweights[i][j])))
                {
                    bestval = fabs(0.5 - arcweights[i][j]);
                    start = i;
                    end = j;
                }
            }
        }
        SCIPfreeBlockMemoryArray( scip, &arcweights[i], modeldata->nC);
    }
    SCIPfreeBlockMemoryArray(scip, &arcweights, modeldata->nC);

    assert(start >= 0 && start <= modeldata->nC - 1);
    assert(end >= 0 && end <= modeldata->nC - 1);

    /* create new nodes */
    SCIP_CALL( SCIPcreateChild(scip, &childprohibit, 0.0, SCIPgetLocalTransEstimate(scip)) );
    SCIP_CALL( SCIPcreateChild(scip, &childenforce, 0.0, SCIPgetLocalTransEstimate(scip)) );

    /* create corresponding constraints */
    SCIP_CALL( SCIPcreateConsArcFlow(scip, &consprohibit, "prohibit", start, end, PROHIBIT, childprohibit, TRUE) );
    SCIP_CALL( SCIPcreateConsArcFlow(scip, &consenforce, "enforce", start, end, ENFORCE, childenforce, TRUE) );

    /* add constraints to nodes */
    SCIP_CALL( SCIPaddConsNode(scip, childprohibit, consprohibit, NULL) );
    SCIP_CALL( SCIPaddConsNode(scip, childenforce, consenforce, NULL) );

    /* release constraints */
    SCIP_CALL( SCIPreleaseCons(scip, &consprohibit) );
    SCIP_CALL( SCIPreleaseCons(scip, &consenforce) );

    *result = SCIP_BRANCHED;
    return SCIP_OKAY;
}

/** creates the arc flow branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleArcFlow(
        SCIP*                 scip                /**< SCIP data structure */
)
{
    SCIP_BRANCHRULEDATA* branchruledata;
    SCIP_BRANCHRULE* branchrule;

    /* create arc flow branching rule data */
    branchruledata = NULL;
    branchrule = NULL;
    /* include branching rule */
    SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC, BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH,
                                          BRANCHRULE_MAXBOUNDDIST, branchruledata) );
    assert(branchrule != NULL);

    SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpArcFlow) );

    return SCIP_OKAY;
}