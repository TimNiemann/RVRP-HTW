/**@file   tools_evaluating.c
 * @brief  functions for evaluation solutions
 * @author Lukas Schürmann University Bonn
 */

#include "tools_evaluating.h"
#include "tools_data.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <scip/scip.h>

/*#define SCIP_DEBUG*/ /* uncomment for some (debugging) output */

#define DISABLE_ALL_LASTSERVICEDATES FALSE /**< (TRUE/FALSE) toggle disabling last service date data (-> more flexible planning) */

static
int getTravelTime(
        model_data*     modeldata,
        int             start,
        int             end,
        int             starttime
)
{
    int traveltime;

    assert(modeldata != NULL);
    assert(0 <= start && start < modeldata->nC);
    assert(0 <= end && end < modeldata->nC);
    if (starttime < NOON_START)
    {
        assert(modeldata->t_travelAM != NULL);
        traveltime = modeldata->t_travelAM[start][end];
    }
    else if (starttime < NOON_END)
    {
        assert(starttime >= NOON_START);
        assert(modeldata->t_travelNoon != NULL);
        traveltime = modeldata->t_travelNoon[start][end];
    }
    else if (starttime < EVENING_START)
    {
        assert(starttime >= NOON_END);
        assert(modeldata->t_travelPM != NULL);
        traveltime = modeldata->t_travelPM[start][end];
    } else {
        assert(starttime >= EVENING_START);
        traveltime = modeldata->t_travel[start][end];
    }
    return traveltime;
}

static
SCIP_RETCODE simulateTour(
        SCIP*               scip,
        model_data*         modelData,
        solutionWindow**    windows,
        int*                simulatedtimes,
        int*                tour,
        int                 length,
        double              devfac,
        int*                totalservice,
        int*                totaltravel,
        int*                totalwaiting,
        int*                numwaits,
        int*                delay,
        int*                numdelays,
        int*                failtime,
        int*                numfails
        ){
    double rnd;
    int i;
    int u, v;
    int service, travel;
    *totalservice = 0;
    *totaltravel = 0;
    *totalwaiting = 0;
    *numwaits = 0;
    *delay = 0;
    *numdelays = 0;
    *failtime = 0;
    *numfails = 0;

    /* depot to first customer */
    u = modelData->nC - 1;
    v = tour[0];
    rnd = devfac * (double) rand() / RAND_MAX;
    travel = rnd * modelData->t_travel_maxDev[u][v];
    simulatedtimes[0] = windows[v]->start_t + travel;
    *totaltravel += modelData->t_travel[u][v] + travel;
    if(simulatedtimes[0] > windows[v]->end_t)
    {
        *delay += simulatedtimes[0] - windows[v]->end_t;
    }
    for(i = 1; i <= length; i++)
    {
        u = tour[i-1];
        if(i < length){
            v = tour[i];
        }else{
            v = modelData->nC - 1;
        }

        rnd = devfac * (double) rand() / RAND_MAX;
        if(RAND_ENCOUNTER && rnd > windows[u]->modelwindow->weigth)
        {
            service = SERVICE_LB;
        }else{
            service = modelData->t_service[u];
            if(RAND_SERVICE)
            {
                service += rnd * modelData->t_service_maxDev[u];
            }
        }
        *totalservice += service;
        travel = getTravelTime(modelData, u, v, simulatedtimes[i - 1] + service);
        if(RAND_TRAVEL)
        {
            rnd = devfac * (double) rand() / RAND_MAX;
            travel += rnd * modelData->t_travel_maxDev[u][v];
        }
        *totaltravel += travel;

        simulatedtimes[i] = simulatedtimes[i-1] + service + travel;
        if(i == length) continue;
        /* if we arrive earlier than actual possible it must be due to an hard time window */
        if(simulatedtimes[i] < windows[v]->start_t)
        {
            *totalwaiting += windows[v]->start_t - simulatedtimes[i];
            (*numwaits)++;
            simulatedtimes[i] = windows[v]->start_t;
        }
        if(simulatedtimes[i] > windows[v]->end_t)
        {
            *delay += simulatedtimes[i] - windows[v]->end_t;
            (*numdelays)++;
        }
        if(simulatedtimes[i] > windows[v]->modelwindow->end_t)
        {
            *failtime += simulatedtimes[i] - windows[v]->modelwindow->end_t;
            (*numfails)++;
        }
    }

    return SCIP_OKAY;
}

static
SCIP_RETCODE gethourminsec(
        int time,
        int* h,
        int* min,
        int* sec
        ){
    *h = time / 3600;
    *min = (time - *h * 3600) / 60;
    *sec = time - *h * 3600 - *min * 60;
    return SCIP_OKAY;
}

static
SCIP_RETCODE printtimes(
        int         totalduration,
        int         service,
        int         travel,
        int         waiting,
        double      numwaits,
        int         failtime,
        double      numfails,
        int         delay,
        int         numcustomers
        ){
    int h, min, sec;
    SCIP_CALL(gethourminsec(totalduration, &h, &min, &sec));
    printf("avg total duration %dh %dmin %dsec\n",h,min,sec);
    SCIP_CALL(gethourminsec(service, &h, &min, &sec));
    printf("avg service time %dh %dmin %dsec (%.2f%%), ", h,min,sec, (double) service * 100 / totalduration);
    SCIP_CALL(gethourminsec(travel, &h, &min, &sec));
    printf("avg travel time %dh %dmin %dsec (%.2f%%), ", h,min,sec, (double) travel * 100 / totalduration);
    SCIP_CALL(gethourminsec(waiting, &h, &min, &sec));
    printf("avg waiting time %dh %dmin %dsec (%.2f%%) - %.2f waiting stops\n",
           h,min,sec, (double) waiting * 100 / totalduration, numwaits);
    SCIP_CALL(gethourminsec(failtime, &h, &min, &sec));
    printf("avg failing time %dh %dmin %dsec (per customer %dmin %dsec) - %.2f failed stops\n",
           h,min,sec, (failtime / numcustomers) / 60, (failtime / numcustomers) - ((int) (failtime / numcustomers)/60)*60, numfails);
    SCIP_CALL(gethourminsec(delay, &h, &min, &sec));
    printf("avg delay %dh %dmin %dsec (per customer %dmin %dsec)\n",h,min,sec,
           (delay / numcustomers) / 60, (delay / numcustomers) - ((int) (delay / numcustomers)/60)*60 );

    return SCIP_OKAY;
}

static
SCIP_RETCODE writeSolution(
        SCIP*            scip,
        char*           outputFile,
        solution_data*  solutionData,
        int             numsims,
        double*         avgdelay,
        double*         avgservicetime,
        double*         avgtraveltime,
        double*         avgwaitingtime,
        double*         avgfailtime,
        double*         avgnumwaits,
        double*         encounterprobs,
        double*         avgnumdelays,
        double*         avgnumfails
        ){
    printf("Write solution in file %s\n",outputFile);

    int total, service, travel, waiting, delay, failtime;
    double numwaits, ep, numdelays, numfails;
    int i;
    int sum;
    FILE *fp;

    service = 0;
    travel = 0;
    waiting = 0;
    delay = 0;
    numwaits = 0.0;
    ep = 0.0;
    numdelays = 0.0;
    failtime = 0;
    numfails = 0.0;
    for(i = 0; i < solutionData->nDays; i++)
    {
        service += (int) avgservicetime[i];
        travel += (int) avgtraveltime[i];
        waiting += (int) avgwaitingtime[i];
        delay += (int) avgdelay[i];
        numwaits += avgnumwaits[i];
        ep += encounterprobs[i] * solutionData->tourLength[i];
        numdelays += avgnumdelays[i];
        failtime += (int) avgfailtime[i];
        numfails += avgnumfails[i];
    }
    total = service + travel + waiting;

    fp = fopen(outputFile, "w");
    if (fp == NULL)
    {
        SCIPwarningMessage(scip, "Can't open outputfile: %s\n", outputFile);
        return SCIP_WRITEERROR;
    }
    fprintf(fp, "{\n");
    fprintf(fp, "\t\"number of days\": %d,\n", solutionData->nDays);
    fprintf(fp, "\t\"number of customers\": %d,\n", solutionData->nC - 1);
    fprintf(fp, "\t\"robustness degree\": %d,\n", solutionData->gamma);
    fprintf(fp, "\t\"window length\": %d,\n", solutionData->windowLength);
    fprintf(fp, "\t\"alpha_one\": %.2f,\n", solutionData->alphas[0]);
    fprintf(fp, "\t\"alpha_two\": %.2f,\n", solutionData->alphas[1]);
    fprintf(fp, "\t\"alpha_three\": %.2f,\n", solutionData->alphas[2]);
    fprintf(fp, "\t\"number of simulations\": %d,\n", numsims);
    fprintf(fp, "\t\"evaluation of days\": [\n");
    for(i = 0; i < solutionData->nDays; i++)
    {
        sum = avgwaitingtime[i] + avgtraveltime[i] + avgservicetime[i];
        fprintf(fp, "\t\t{\n");
        fprintf(fp, "\t\t\t\"day\": %d,\n", i);
        fprintf(fp, "\t\t\t\"customers on this day\": %d,\n", solutionData->tourLength[i]);
        fprintf(fp, "\t\t\t\"average total duration\": %d,\n", sum);

        fprintf(fp, "\t\t\t\"average service time\":{\n");
        fprintf(fp, "\t\t\t\t\"total\": %d,\n",(int) avgservicetime[i]);
        fprintf(fp, "\t\t\t\t\"proportion\": %.2f\n",(double) avgservicetime[i]/sum);
        fprintf(fp, "\t\t\t},\n");

        fprintf(fp, "\t\t\t\"average travel time\":{\n");
        fprintf(fp, "\t\t\t\t\"total\": %d,\n",(int) avgtraveltime[i]);
        fprintf(fp, "\t\t\t\t\"proportion\": %.2f\n",(double) avgtraveltime[i]/sum);
        fprintf(fp, "\t\t\t},\n");

        fprintf(fp, "\t\t\t\"average waiting time\":{\n");
        fprintf(fp, "\t\t\t\t\"total\": %d,\n",(int) avgwaitingtime[i]);
        fprintf(fp, "\t\t\t\t\"proportion\": %.2f,\n",(double) avgwaitingtime[i]/sum);
        fprintf(fp, "\t\t\t\t\"avg num of waits\": %.2f\n", avgnumwaits[i]);
        fprintf(fp, "\t\t\t},\n");

        fprintf(fp, "\t\t\t\"average delay time\":{\n");
        fprintf(fp, "\t\t\t\t\"total\": %d,\n",(int) avgdelay[i]);
        fprintf(fp, "\t\t\t\t\"per customer\": %.2f,\n",(double) avgdelay[i]/solutionData->tourLength[i]);
        fprintf(fp, "\t\t\t\t\"avg num of delays\": %.2f\n", avgnumdelays[i]);
        fprintf(fp, "\t\t\t},\n");

        fprintf(fp, "\t\t\t\"average fail time\": %d,\n", (int) avgfailtime[i]);
        fprintf(fp, "\t\t\t\"average num of fails\": %.2f,\n", avgnumfails[i]);
        fprintf(fp, "\t\t\t\"encounter probability\": %.2f\n", encounterprobs[i]);
        fprintf(fp, "\t\t}");
        if(i != solutionData->nDays - 1)
        {
            fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\t],\n");
    fprintf(fp, "\t\"total evaluation\": {\n");

    fprintf(fp, "\t\t\"average total duration\": %d,\n", total);

    fprintf(fp, "\t\t\"average service time\":{\n");
    fprintf(fp, "\t\t\t\"total\": %d,\n",service);
    fprintf(fp, "\t\t\t\"proportion\": %.2f\n",(double) service/total);
    fprintf(fp, "\t\t},\n");

    fprintf(fp, "\t\t\"average travel time\":{\n");
    fprintf(fp, "\t\t\t\"total\": %d,\n",travel);
    fprintf(fp, "\t\t\t\"proportion\": %.2f\n",(double) travel/total);
    fprintf(fp, "\t\t},\n");

    fprintf(fp, "\t\t\"average waiting time\":{\n");
    fprintf(fp, "\t\t\t\"total\": %d,\n",waiting);
    fprintf(fp, "\t\t\t\"proportion\": %.2f,\n",(double) waiting/total);
    fprintf(fp, "\t\t\t\"avg num of waits\": %.2f\n", numwaits);
    fprintf(fp, "\t\t},\n");

    fprintf(fp, "\t\t\"average delay time\":{\n");
    fprintf(fp, "\t\t\t\"total\": %d,\n",delay);
    fprintf(fp, "\t\t\t\"per customer\": %.2f,\n",(double) delay / (solutionData->nC - 1));
    fprintf(fp, "\t\t\t\"avg num of delays\": %.2f\n", numdelays);
    fprintf(fp, "\t\t},\n");

    fprintf(fp, "\t\t\t\"average fail time\": %d,\n", failtime);
    fprintf(fp, "\t\t\"average num of fails\": %.2f,\n", numfails);
    fprintf(fp, "\t\t\"encounter probability\": %.2f\n", ep / (solutionData->nC - 1));

    fprintf(fp, "\t}\n");

    fprintf(fp, "}");
    fclose(fp);
    return SCIP_OKAY;
}

SCIP_RETCODE evaluateSolution(
        SCIP*           scip,
        model_data*     modelData,
        solution_data*  solutionData,
        char*           outputFile,
        double          devfac,
        int             numsims,
        SCIP_Bool       verbose
        ){
    FILE* fp;

    int i, node;
    int day;
    int* simulatedtimes;
    int* tour;
    int length;

    int servicetime, traveltime, waitingtime, failtime;
    int numwaits;
    int numdelays;
    int numfails;

    double *avgsericetime, *avgtraveltime, *avgwaitingtime, *avgfailtime;
    int sum;
    double *avgnumwaits;
    double *avgnumdelays;
    double *avgnumfails;

    int delay;
    int sim;

    solutionWindow* solutionwindow;

    double* avgdelay;
    double* encounterprobs;
    double totalencounterprob;

    int totalavgduration = 0;
    int totalavgservice = 0;
    int totalavgtravel = 0;
    int totalavgwaiting = 0;
    double totalavgnumwaits = 0;
    int totalavgdelay = 0;
    double totalavgnumdelays = 0;
    int totalavgfailtime = 0;
    double totalavgnumfails = 0;

    srand(RNDSEED);

    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &simulatedtimes, modelData->nC));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &encounterprobs, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgwaitingtime, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgfailtime, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgtraveltime, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgsericetime, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgnumwaits, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgnumfails, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgdelay, modelData->nDays));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip, &avgnumdelays, modelData->nDays));
    if(verbose)
    {
        printf("%d Customers, %d days\n", modelData->nC - 1, modelData->nDays);
        printf("%d simulations\n", numsims);
        printf("Deviation factor: %.2f\n", devfac);
    }
    for(day = 0; day < modelData->nDays; day++)
    {
        avgsericetime[day] = 0.0;
        avgtraveltime[day] = 0.0;
        avgwaitingtime[day] = 0.0;
        avgfailtime[day] = 0.0;
        avgnumwaits[day] = 0.0;
        avgnumfails[day] = 0.0;
        encounterprobs[day] = 0.0;
        avgdelay[day] = 0.0;
        avgnumdelays[day] = 0.0;
        tour = solutionData->tour[day];
        length = solutionData->tourLength[day];
        if(length == 0) continue;

        for(i = 0; i < length; i++)
        {
            node = tour[i];
            solutionwindow = solutionData->timeWindows[node];
            encounterprobs[day] += solutionwindow->modelwindow->weigth / length;
        }
        if(verbose)
        {
            printf("---------------- Day %d ----------------\n", day);
            printf("%d customers, avg encounter probability %.2f%%\n", length, encounterprobs[day] * 100);
        }
        /** Simulate #numsims tours for the current day */
        for(sim = 0; sim < numsims; sim++)
        {
            SCIP_CALL( simulateTour(scip, modelData, solutionData->timeWindows, simulatedtimes, tour, length, devfac,
                                    &servicetime, &traveltime, &waitingtime, &numwaits, &delay, &numdelays, &failtime, &numfails));
            avgsericetime[day] += (double) servicetime / numsims;
            avgtraveltime[day] += (double) traveltime / numsims;
            avgwaitingtime[day] += (double) waitingtime / numsims;
            avgfailtime[day] += (double) failtime / numsims;
            avgnumwaits[day] += (double) numwaits / numsims;
            avgnumfails[day] += (double) numfails / numsims;
            avgdelay[day] += (double) delay / numsims;
            avgnumdelays[day] += (double) numdelays / numsims;
        }
        sum = avgwaitingtime[day] + avgtraveltime[day] + avgsericetime[day];
        if(verbose) printtimes(sum, (int) avgsericetime[day], (int) avgtraveltime[day], (int) avgwaitingtime[day],
                               avgnumwaits[day], (int) avgfailtime[day], avgnumfails[day], (int) avgdelay[day], length);
        totalavgduration += sum;
        totalavgservice += avgsericetime[day];
        totalavgtravel += avgtraveltime[day];
        totalavgwaiting += avgwaitingtime[day];
        totalavgnumwaits += avgnumwaits[day];
        totalavgfailtime += avgfailtime[day];
        totalavgnumfails += avgnumfails[day];
        totalavgdelay += avgdelay[day];
        totalavgnumdelays += avgnumdelays[day];
    }
    if(verbose) printf("---------------- TOTAL ----------------\n");

    totalencounterprob = 0.0;
    for(day = 0; day < modelData->nDays; day++)
    {
        totalencounterprob += encounterprobs[day] * solutionData->tourLength[day];
    }
    fp = fopen("evaluation.out", "a");
    fprintf(fp, "%s & %d & %d & %.2f & %.2f & %.2f & %d & %.2f & %d & %.2f & %d & %.2f & %d & %d & %.2f & %.2f & %d & %.2f & %.2f\n",
           outputFile,
           solutionData->gamma, solutionData->windowLength, solutionData->alphas[0], solutionData->alphas[1], solutionData->alphas[2], 
           totalavgservice, (double) 100* totalavgservice/totalavgduration,
           totalavgtravel, (double) 100* totalavgtravel/totalavgduration,
           totalavgwaiting, (double) 100* totalavgwaiting/totalavgduration,
           totalavgduration,
           totalavgdelay, (double) totalavgdelay/(modelData->nC-1), totalavgnumdelays,
           totalavgfailtime, totalavgnumfails,
           totalencounterprob * 100 / (modelData->nC -1));
    if(verbose)
    {
        printf("total avg encounter probability %.2f%%\n",totalencounterprob * 100 / modelData->nC);
        printtimes(totalavgduration, totalavgservice,totalavgtravel,totalavgwaiting,totalavgnumwaits,totalavgfailtime,totalavgnumfails,totalavgdelay,modelData->nC - 1);
    }
    if(outputFile != NULL) {
        SCIP_CALL(writeSolution(scip, outputFile, solutionData, numsims,
                                avgdelay, avgsericetime, avgtraveltime, avgwaitingtime, avgfailtime, avgnumwaits, encounterprobs, avgnumdelays, avgnumfails));
    }
    double maxep = 0.0;
    modelWindow* window;
    double currentmax;
    for(i = 0; i < modelData->nC - 1; i++)
    {
        currentmax = 0.0;
        window = modelData->timeWindows[i];
        while(window != NULL)
        {
            if(currentmax < window->weigth)
            {
                currentmax = window->weigth;
            }
            window = window->next;
        }
        maxep += currentmax / modelData->nC;
    }
    if(verbose) printf("max ep = %.2f%%\n",100 * maxep);

    SCIPfreeBlockMemoryArray(scip, &avgdelay, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgsericetime, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgtraveltime, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgwaitingtime, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgnumwaits, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgfailtime, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &avgnumfails, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &encounterprobs, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &simulatedtimes, modelData->nC);
    return SCIP_OKAY;
}
