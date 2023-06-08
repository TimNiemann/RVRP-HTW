/**@file   postprocessing_vrp.c
 * @brief  solution structure and functions to create an output file
 * @author Tim Niemann, TU Braunschweig
 */

#include <stdio.h>
#include "scip/scip.h"
#include "postprocessing_vrp.h"
#include "tools_data.h"
#include "tools_vrp.h"
#include "vardata_vrp.h"
#include "probdata_vrp.h"
#include "tools_evaluating.h"

extern
SCIP_RETCODE createSolutionWindow(
    SCIP*               scip,
    solutionWindow**    solutionwindow,
    int                 day,
    int                 start,
    int                 end,
    int                 latestArrival,
    modelWindow*        window
    )
{
    assert(solutionwindow != NULL);
    assert(window != NULL);
    SCIP_CALL( SCIPallocMemory(scip, solutionwindow) );
    (*solutionwindow)->day = day;
    (*solutionwindow)->start_t = start;
    (*solutionwindow)->end_t = end;
    (*solutionwindow)->latest_arrival = latestArrival;
    SCIP_CALL( SCIPduplicateMemory(scip, &(*solutionwindow)->modelwindow, window));
    assert((*solutionwindow)->modelwindow != NULL);
    return SCIP_OKAY;
}

extern
SCIP_RETCODE freeSolutionWindow(
    SCIP*               scip,
    solutionWindow*     solutionwindow
    )
{
    if (solutionwindow != NULL)
    {
        SCIPfreeMemory(scip, &(solutionwindow->modelwindow));
        SCIPfreeMemory(scip, &solutionwindow);
    }
    return SCIP_OKAY;
}

extern
SCIP_RETCODE freeSolutionWindowArray(
    SCIP*               scip,
    solutionWindow**    solutionwindow,
    int                 size
    )
{
    int i;
    if (solutionwindow != NULL)
    {
        for (i = 0; i < size; i++)
        {
            if(solutionwindow[i] == NULL) printf("i %d, size %d\n",i,size);
            assert(solutionwindow[i] != NULL);
            SCIP_CALL( freeSolutionWindow(scip, solutionwindow[i]) );
        }
        SCIPfreeMemoryArray(scip, &solutionwindow);
    }
    return SCIP_OKAY;
}

static
SCIP_RETCODE createSolutionData(
    SCIP*               scip,
    solution_data**     solutionData,
    model_data*         modeldata,
    solutionWindow**    timeWindows,
    SCIP_Real*          obj,
    int**               tour,
    int*                tourLength,
    int*                tourDuration
    )
{
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    int i;
    assert(scip != NULL);
    assert(solutionData != NULL);
    assert(*solutionData == NULL);
    assert(modeldata != NULL);

    SCIP_CALL (SCIPallocMemory(scip, solutionData) );
    (*solutionData)->nC = modeldata->nC;
    (*solutionData)->nDays = modeldata->nDays;
    (*solutionData)->gamma = modeldata->maxDelayEvents;
    (*solutionData)->windowLength = probdata->delayTolerance;
    (*solutionData)->objCoefficient = probdata->weightObjective;
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->timeWindows, timeWindows, modeldata->nC - 1) );
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->obj, obj, modeldata->nDays) );
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->tourLength, tourLength, modeldata->nDays) );
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->tourDuration, tourDuration, modeldata->nDays) );
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->customerIDs, modeldata->customerIDs, modeldata->nC) );
    SCIP_CALL (SCIPduplicateMemoryArray(scip, &(*solutionData)->alphas, probdata->alphas, 3) );
    SCIP_CALL (SCIPallocMemoryArray(scip, &((*solutionData)->tour), modeldata->nDays) );
    for (i = 0; i < modeldata->nDays; i++)
    {
        SCIP_CALL( SCIPduplicateMemoryArray(scip, &(*solutionData)->tour[i], tour[i], tourLength[i]) );
    }
    for (i = 0; i < modeldata->nC - 1; i++)
    {
        SCIP_CALL( createSolutionWindow(scip, &(*solutionData)->timeWindows[i], timeWindows[i]->day, timeWindows[i]->start_t, timeWindows[i]->end_t, timeWindows[i]->latest_arrival, timeWindows[i]->modelwindow) );
    }
    return SCIP_OKAY;
}

static
SCIP_RETCODE extractSolution(
    SCIP*               scip,
    solution_data**     solutionData,
    model_data*         modeldata
    )
{
    SCIP_STAGE stage;
    SCIP_SOL* sol = NULL;
    SCIP_VAR** vars = NULL;
    SCIP_VARDATA** vardatas = NULL;
    SCIP_PROBDATA* probdata = SCIPgetProbData(scip);
    double* vals = NULL;
    int nvars, nbinvars;
    int i, j = 0;

    solutionWindow** timeWindows = NULL;
    double* obj = NULL;
    int** tour;
    int* tourLength = NULL;
    int* tourDuration = NULL;

    assert(scip != NULL);

    stage = SCIPgetStage(scip);
    if (stage != SCIP_STAGE_SOLVING && stage != SCIP_STAGE_SOLVED)
    {
        SCIPwarningMessage(scip, "Method can not be called in this stage: %d\n", stage);
        return SCIP_INVALIDCALL;
    }
    if (SCIPgetNSols(scip) == 0)
    {
        SCIPwarningMessage(scip, "SCIP can't find a feasible solution. No output file was created.\n");
        return SCIP_FILECREATEERROR;
    }

    sol = SCIPgetBestSol(scip);
    assert(sol != NULL);

    SCIP_CALL( SCIPgetSolVarsData(scip, sol, &vars, &nvars, &nbinvars, NULL, NULL, NULL) );
    assert(nvars == nbinvars);

    SCIP_CALL( SCIPallocBufferArray(scip, &vals, nvars) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &vardatas, modeldata->nDays) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &obj, modeldata->nDays) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &tourDuration, modeldata->nDays) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &tourLength, modeldata->nDays) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &tour, modeldata->nDays) );
    SCIP_CALL( SCIPallocMemoryArray(scip, &timeWindows, modeldata->nC - 1) );

    for (i = 0; i < modeldata->nC -1; i++)
    {
        timeWindows[i] = NULL;
    }

    SCIP_CALL (SCIPgetSolVals(scip, sol, nvars, vars, vals) );
    for (i = 0; i < nvars; i++)
    {
        if (SCIPisEQ(scip, vals[i], 1.0))
        {
            int k;
            j = SCIPvarGetData(vars[i])->day;
            vardatas[j] = SCIPvarGetData(vars[i]);
            obj[j] = SCIPvarGetObj(vars[i]);
            tourDuration[j] = vardatas[j]->tourduration;
            tourLength[j] = vardatas[j]->tourlength;
            SCIP_CALL( SCIPallocMemoryArray(scip, &tour[j], tourLength[j]) );
            for (k = 0; k < tourLength[j]; k++)
            {
                int start, end, latestArrival;
                modelWindow* window;
                tour[j][k] = vardatas[j]->customertour[k];
                assert(!nodeIsDepot(modeldata, tour[j][k]));
                /* create time windows for delivery */
                start = vardatas[j]->solutionwindows[k]->start_t;
                end = vardatas[j]->solutionwindows[k]->end_t;
                latestArrival = vardatas[j]->solutionwindows[k]->latest_arrival;
                window = vardatas[j]->solutionwindows[k]->modelwindow;
                SCIP_CALL( createSolutionWindow(scip, &timeWindows[tour[j][k]], vardatas[j]->day, start, end, latestArrival, window) );
            }
        }
    }

    /* create empty entries for optional customers because of price collecting */
    printf("\n");
    for (i = 0; i < modeldata->nC -1; i++)
    {
        if (probdata->useOptionals == TRUE && probdata->optionalCustomers[i] == TRUE && timeWindows[i] == NULL)
        {
            modelWindow* window = NULL;
            SCIP_CALL( SCIPallocMemory(scip, &window) );
            window->day = -1;
            window->end_t = -1;
            window->start_t = -1;
            window->weigth = 0;
            window->next = NULL;
            SCIP_CALL( createSolutionWindow(scip, &timeWindows[i], -1, -1, -1, -1, window) );
            printf(" * Optional customer %d is not served!\n", i);
        }
    }
    
    SCIP_CALL( createSolutionData(scip, solutionData, modeldata, timeWindows, obj, tour, tourLength, tourDuration) );

    /* free memory */
    SCIP_CALL( freeSolutionWindowArray(scip, timeWindows, modeldata->nC - 1) );
    for (i = 0; i < modeldata->nDays; i++)
    {
        SCIPfreeMemoryArray(scip, &tour[i]);
    }
    SCIPfreeMemoryArray(scip, &tour);
    SCIPfreeMemoryArray(scip, &tourLength);
    SCIPfreeMemoryArray(scip, &tourDuration);
    SCIPfreeMemoryArray(scip, &obj);
    SCIPfreeMemoryArray(scip, &vardatas);
    SCIPfreeBufferArray(scip, &vals);    
    return SCIP_OKAY;
}

static
SCIP_RETCODE freeSolutionData(
    SCIP*               scip,
    solution_data*      solutionData
    )
{
    int i;
    assert(solutionData != NULL);
    
    SCIPfreeMemoryArray(scip, &solutionData->obj);
    SCIPfreeMemoryArray(scip, &solutionData->tourDuration);
    SCIPfreeMemoryArray(scip, &solutionData->customerIDs);
    for (i = 0; i < solutionData->nDays; i++)
    {
        SCIPfreeMemoryArray(scip, &solutionData->tour[i]);
    }
    SCIPfreeMemoryArray(scip, &solutionData->tourLength);
    SCIPfreeMemoryArray(scip, &solutionData->tour);
    SCIPfreeMemoryArray(scip, &solutionData->alphas);
    SCIP_CALL( freeSolutionWindowArray(scip, solutionData->timeWindows, solutionData->nC - 1) );
    SCIPfreeMemory(scip, &solutionData);

    return SCIP_OKAY;
}


/**
 * @return average Delay of the soltion in seconds */
static
void printSolutionStatistics(
    solution_data*      solution
    )
{
    int averageDelay = 0.0;
    int averageDuration = 0.0;
    int maxDelay = 0;
    int maxDuration = 0;
    int usedDays = 0;
    int tmp;
    int i;
    assert(solution != NULL);

    for (i = 0; i < solution->nC - 1; i++)
    {
        /* delay */
        tmp = solution->timeWindows[i]->latest_arrival - solution->timeWindows[i]->start_t;
        averageDelay += tmp;
        if (tmp > maxDelay)
        {
            maxDelay = tmp;
        }
    }
    for (i = 0; i < solution->nDays; i++)
    {
        /* duration */
        tmp = solution->tourDuration[i];
        averageDuration += tmp;
        if (tmp > maxDuration)
        {
            maxDuration = tmp;
        }
        if (tmp > 0)
        {
            usedDays++;
        }
    }
    averageDelay /= solution->nC - 1;
    averageDuration /= usedDays;

    /* convert duration to minutes */
    averageDuration /= 60;
    maxDuration /= 60;

    printf("Delay of Solution: Average: %d:%02d min, Maximum: %d:%02d min.\n", averageDelay/60, averageDelay%60, maxDelay/60, maxDelay%60);
    printf("Duration of Solution: Average: %d:%02d h, Maximum: %d:%02d h.\n", averageDuration/60, averageDuration%60, maxDuration/60, maxDuration%60);
}

extern
SCIP_RETCODE createOutputFile(
    SCIP*               scip,
    model_data*         modeldata,
    char*               outputFile
    )
{
    solution_data* solutionData = NULL;
    assert(scip != NULL);
    assert(modeldata != NULL);
    assert(outputFile != NULL);

    SCIP_CALL( extractSolution(scip, &solutionData, modeldata) ); 
    assert(solutionData != NULL);
    SCIP_CALL( writeSolutionJSON(scip, solutionData, outputFile) );
    printSolutionStatistics(solutionData);

    SCIP_CALL( freeSolutionData(scip, solutionData) );
    return SCIP_OKAY;
}