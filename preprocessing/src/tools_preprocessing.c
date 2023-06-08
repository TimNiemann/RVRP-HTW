/**@file   tools_preprocessing.c
 * @brief  functions for handling of the preprocessing of the instances
 * @author Lukas Schürmann University Bonn
 */

#include "tools_preprocessing.h"
#include "tools_data.h"
#include <string.h>
#include <math.h>
#include <scip/scip.h>

/*#define SCIP_DEBUG*/ /* uncomment for some (debugging) output */

#define DISABLE_ALL_LASTSERVICEDATES FALSE /**< (TRUE/FALSE) toggle disabling last service date data (-> more flexible planning) */


/** read command line arguments */
SCIP_RETCODE readArguments(
        int          argc,                /**< number of shell parameters */
        char**       argv,                /**< array with shell parameters */
        char**       dataFile,            /**< path of data file (e.g., ../inst_bayern.dat) */
        char**       outputFile,          /**< optional: output instance file (incl. path; default: inst_<dataSetName>.dat) */
        int          *maxWaitingTime,     /**< optional: maximum waiting time between two customers */
        int          *maxDays,            /**< optional: maximum amount of days to be considered */
        int          *maxTime,            /**< optional: ignore everything after this point of time */
        int          *reducednC,          /**< optional: maximum amount of customers to be considered */
        double       *windowPercentage    /**< optional: percentage of time windows to be removed */
        )
{
    char usage[SCIP_MAXSTRLEN];
    char strbuffer[SCIP_MAXSTRLEN];
    char* locstr;
    int status;

    assert( argc >= 1 );
    assert( argv != NULL );
    assert( dataFile != NULL );
    assert( outputFile != NULL );

    /* init usage text */
    status = snprintf(usage, SCIP_MAXSTRLEN - 1,
                      "usage: %s <path of file> \
                      [-o <output instance file (optional; default used if no name is provided or if it begins with \"-\")>] \
                      [-h <maximum waiting time between two customers in seconds (optional; default : 3600) >] \
                      [-t <upper bound for daily working time (mandatory; if = 0: create full flexibility for each customer)>] \
                      [-d <number of considered days (mandatory)>] \
                      [-r <number of considered customers (mandatory)>] \
                      [-w <percentage of time windows to be removed (mandatory, in [0,1])>]", argv[0]);
    assert( 0 <= status && status < SCIP_MAXSTRLEN );

    /* init arguments */
    *dataFile           = NULL;
    *outputFile          = NULL;

    /* first mandatory argument: path to base data folder */
    *dataFile = argv[1];
    if ( *dataFile == NULL )
    {
        fprintf(stderr, "No path of data supplied.\n");
        fprintf(stderr, "%s\n", usage);
        return SCIP_ERROR;
    }

    /* check for optional arguments */
    for (int i = 2; i < argc; i++)
    {
        if ( ! strcmp(argv[i], "-o") )
        {
            locstr = (char*) malloc((int)(strlen(argv[1])-3)*sizeof(char));
            snprintf(locstr, (int) strlen(argv[1])-3,"%s", argv[1]);
            assert( *outputFile == NULL );
            *outputFile = (char *) malloc ( SCIP_MAXSTRLEN * sizeof(char) +1 );

            if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
            {
                sprintf(strbuffer, "%s_preprocessed.dat", locstr);
                strcpy(*outputFile,strbuffer);
            }
            else
            {
                i++;
                strcpy(*outputFile,argv[i]);
                assert( i < argc );
            }
            assert( *outputFile != NULL );
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-h"))
        {
            if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
            {
                *maxWaitingTime = 3600;
            }
            else
            {
                i++;
                locstr = (char*) malloc((int)(strlen(argv[i]))*sizeof(char));
                snprintf(locstr, (int) strlen(argv[i]),"%s", argv[i]);
                *maxWaitingTime = atoi(locstr);
                free(locstr);
            }
        }
        if ( ! strcmp(argv[i], "-t") )
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing Max. working time. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *maxTime = atoi(locstr);
            if(*maxTime != 0 && (*maxTime < DEFAULT_SHIFT_START || *maxTime > DEFAULT_SHIFT_END ))
            {
                fprintf(stderr, "Invalid Max. working time -> Choose from interval (%d, %d) or set to 0. ",DEFAULT_SHIFT_START,DEFAULT_SHIFT_END);
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-d"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing number of days. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *maxDays = atoi(locstr);
            if(*maxDays <= 0 )
            {
                fprintf(stderr, "Invalid number of days -> needs to be larger than 0. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-r"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing reduced number of customer. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *reducednC = atoi(locstr);
            if(*reducednC <= 0 )
            {
                fprintf(stderr, "Invalid number of customers -> needs to be larger than 0. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-w"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing percentage of windows to be removed. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *windowPercentage = atof(locstr);
            if(*windowPercentage < 0.0 || *windowPercentage > 1.0)
            {
                fprintf(stderr, "Invalid percentage of windows to be removed -> needs to be in [0,1]. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }

    }
    return SCIP_OKAY;
}

/** creates a time window from shift start to shift end for each day for each customer */
SCIP_RETCODE fullFlexibility(
        SCIP*           scip,
        model_data*     modelData
){
    int i, j, k;
    modelWindow* window;
    modelWindow* lastwindow;

    assert(modelData->timeWindows != NULL);
    for(i = 0; i < modelData->nC - 1; i++)
    {
        if(modelData->timeWindows[i] != NULL)
        {
            k = 0;
            window = modelData->timeWindows[i];
            while(window != NULL)
            {
                window->start_t = DEFAULT_SHIFT_START;
                window->end_t = DEFAULT_SHIFT_END;
                window->day = k++;
                lastwindow = window;
                window = window->next;
            }
            for(j = k; j < modelData->nDays; j++)
            {
                SCIP_CALL( SCIPallocBlockMemory(scip, &window) );
                window->start_t = DEFAULT_SHIFT_START;
                window->end_t = DEFAULT_SHIFT_END;
                window->day = j;
                window->next = NULL;

                lastwindow->next = window;
                lastwindow = window;
                window = window->next;
            }
        }
    }
    return SCIP_OKAY;
}

/** removes a given percentage of all time windows, randomly chosen.
 *  The chance of removal is applied for every window independently, which creates a binomial distribution.
 */
SCIP_RETCODE removeTimeWindowsBinomial(
        SCIP* scip,
        model_data* modelData,
        double percentage
        )
{
    int i;
    modelWindow* tw = NULL;
    modelWindow* lasttw = NULL;
    assert(0 <= percentage && percentage <= 1);

    for(i=0; i<modelData->nC-1; i++)/* check for each customer */
    {
        if(modelData->timeWindows != NULL)
        {
            lasttw = NULL;
            tw = modelData->timeWindows[i];
            while (tw != NULL) /* check for each of the customers time windows */
            {
                double rnd = (double) rand() / RAND_MAX;
                if(rnd <= percentage) /* if random number is lower than the given percentage, we delete the window */
                {
                    if(lasttw == NULL)
                    {
                        modelData->timeWindows[i] = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = modelData->timeWindows[i];
                    }
                    else
                    {
                        lasttw->next = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = lasttw->next;
                    }
                }
                else
                {
                    lasttw = tw;
                    tw = tw->next;
                }
            }
        }
    }

    return SCIP_OKAY;
}

/** removes a given percentage of all time windows, randomly chosen.
 *  The chance of removal is applied for all windows of a customer together.
 */
SCIP_RETCODE removeTimeWindowsEvenly(
        SCIP* scip,
        model_data* modelData,
        double percentage
        )
{
    int i;
    modelWindow* tw = NULL;
    modelWindow* lasttw = NULL;
    assert(0 <= percentage && percentage <= 1);
    if (modelData->timeWindows == NULL)
    {
        return SCIP_OKAY;
    }

    for(i=0; i<modelData->nC-1; i++)/* check for each customer */
    {
        int nWindows = 0;
        int nDeleted = 0;
        int nRemaining = 0;
        lasttw = NULL;
        tw = modelData->timeWindows[i];
        while(tw != NULL)
        {
            tw = tw->next;
            nWindows++;
        }
        nRemaining = nWindows;

        tw = modelData->timeWindows[i];
        while (tw != NULL && nDeleted < nWindows * percentage) /* check for each of the customers time windows */
        {
            double rnd = (double) rand() / RAND_MAX;
            if(rnd <= percentage || nRemaining - nDeleted <= nWindows * percentage) /* if random number is lower than the given percentage, we delete the window */
            {
                nDeleted++;
                if(lasttw == NULL)
                {
                    modelData->timeWindows[i] = tw->next;
                    SCIPfreeBlockMemory(scip, &tw);
                    tw = modelData->timeWindows[i];
                }
                else
                {
                    lasttw->next = tw->next;
                    SCIPfreeBlockMemory(scip, &tw);
                    tw = lasttw->next;
                }
            }
            else
            {
                lasttw = tw;
                tw = tw->next;
            }
            nRemaining--;
        }
    }

    return SCIP_OKAY;
}

/** deletes/cuts every time window that starts/ends after value maxTime */
SCIP_RETCODE cutTimeWindows(
        SCIP* scip,
        model_data* modelData,
        int         maxTime
        )
{
    int i;
    modelWindow* tw = NULL;
    modelWindow* lasttw = NULL;

    for(i=0; i<modelData->nC-1; i++)/* check for each customer */
    {
        if(modelData->timeWindows != NULL)
        {
            lasttw = NULL;
            tw = modelData->timeWindows[i];
            while (tw != NULL) /* check for each of the customers time windows */
            {
                if(tw->start_t >= maxTime) /* if the time window lies entirely behind maxTime we delete it */
                {
                    if(lasttw == NULL)
                    {
                        modelData->timeWindows[i] = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = modelData->timeWindows[i];
                    }
                    else
                    {
                        lasttw->next = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = lasttw->next;
                    }
                }
                else if(tw->end_t >= maxTime) /* else if the end time lies behind maxTime we cut it*/
                {
                    tw->end_t = maxTime;
                    lasttw = tw;
                    tw = tw->next;
                }
                else
                {
                    lasttw = tw;
                    tw = tw->next;
                }
            }
        }
    }

    return SCIP_OKAY;
}

/** deletes time windows of days after maxDay */
SCIP_RETCODE deleteDays(
        SCIP* scip,
        model_data* modelData,
        int         maxDays
        )
{
    int i;
    modelWindow* tw = NULL;
    modelWindow* lasttw = NULL;

    for(i=0; i<modelData->nC-1; i++)/* check for each customer */
    {
        if(modelData->timeWindows != NULL)
        {
            lasttw = NULL;
            tw = modelData->timeWindows[i];
            while (tw != NULL) /* check for each of the customers time windows */
            {
                if(tw->day >= maxDays) /* if the time window takes place on a day after maxDays */
                {
                    if(lasttw == NULL)
                    {
                        modelData->timeWindows[i] = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = modelData->timeWindows[i];
                    }
                    else
                    {
                        lasttw->next = tw->next;
                        SCIPfreeBlockMemory(scip, &tw);
                        tw = lasttw->next;
                    }
                }
                else
                {
                    assert(tw != NULL);
                    lasttw = tw;
                    tw = tw->next;
                }
            }
        }
    }

    return SCIP_OKAY;
}

/** merges all pairs of successive time windows */
SCIP_RETCODE mergeTimeWindows(
        SCIP* scip,
        model_data* modelData
){
    int i;
    modelWindow* window;
    modelWindow* secondwindow;
    modelWindow* prevwindow;

    for(i = 0; i < modelData->nC - 1; i++)
    {
        window = modelData->timeWindows[i];
        while(window != NULL)
        {
            secondwindow = modelData->timeWindows[i];
            prevwindow = NULL;
            while (secondwindow != NULL)
            {
                if(secondwindow->day == window->day && secondwindow != window)
                {
                    if(window->start_t <= secondwindow->start_t)
                    {
                        if(window->end_t >= secondwindow->start_t)
                        {
                            if(window->end_t <= secondwindow->end_t)
                            {
                                window->end_t = secondwindow->end_t;
                            }
                            window->weigth = (window->weigth + secondwindow->weigth) / 2;
                            if(prevwindow != NULL)
                            {
                                prevwindow->next = secondwindow->next;
                            }else{
                                modelData->timeWindows[i] = secondwindow->next;
                            }
                            SCIPfreeBlockMemory(scip, &secondwindow);
                            secondwindow = prevwindow;
                        }
                    }else{
                        if(secondwindow->end_t >= window->start_t)
                        {
                            window->start_t = secondwindow->start_t;
                            if(window->end_t < secondwindow->end_t)
                            {
                                window->end_t = secondwindow->end_t;
                            }
                            window->weigth = (window->weigth + secondwindow->weigth) / 2;
                            if(prevwindow != NULL)
                            {
                                prevwindow->next = secondwindow->next;
                            }else{
                                modelData->timeWindows[i] = secondwindow->next;
                            }
                            SCIPfreeBlockMemory(scip, &secondwindow);
                            secondwindow = prevwindow;
                        }
                    }
                }
                if(secondwindow != NULL)
                {
                    prevwindow = secondwindow;
                    secondwindow = secondwindow->next;
                }else{
                    secondwindow = modelData->timeWindows[i];
                }
            }
            window = window->next;
        }
    }

    return SCIP_OKAY;
}

/** set up the neighborhoods with help of the time windows */
SCIP_RETCODE processNeighbors(
        SCIP* scip,
        model_data* modelData,
        int deletedCustomer,
        int maxWaitingTime
        )
{
    int i,j,k,l;
    int day;
    int t1;
    int sortvalue;
    int** values;
    int* numneighbors;
    SCIP_Bool** isneighbor;
    SCIP_Bool* repeat;
    modelWindow* tw1 = NULL;
    modelWindow* tw2 = NULL;
    neighbor* n = NULL;
    neighbor* nlast = NULL;

    /* free data if already allocated - for heuristic edit */
    if(modelData->neighbors != NULL)
    {
        for(i = 0; i < modelData->nC+deletedCustomer; i++)
        {
            if(modelData->neighbors[i] != NULL)
            {
                for(j = 0; j < modelData->nDays; j++)
                {
                    while(modelData->neighbors[i][j] != NULL)
                    {
                        n = modelData->neighbors[i][j];
                        modelData->neighbors[i][j] = (modelData->neighbors[i][j])->next;
                        SCIPfreeBlockMemory(scip, &n);
                    }
                }
                SCIPfreeBlockMemoryArray(scip, &(modelData->neighbors[i]), modelData->nDays);
            }
        }
        SCIPfreeBlockMemoryArray(scip, &(modelData->neighbors), modelData->nC+deletedCustomer);
    }
    /* init modeldata->neighbors */
    assert(modelData->neighbors == NULL);
    if(modelData->neighbors == NULL)
    {
        SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->neighbors),modelData->nC+deletedCustomer));
        for(i=0; i<modelData->nC+deletedCustomer; i++)
        {
            SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->neighbors[i]),modelData->nDays));
            for(j=0; j<modelData->nDays; j++)
            {
                modelData->neighbors[i][j] = NULL;
            }
        }
    }
    /* alloc for first neighbor */
    SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
    SCIPallocBlockMemoryArray(scip, &values, modelData->nDays);
    SCIPallocBlockMemoryArray(scip, &isneighbor, modelData->nDays);
    SCIPallocBlockMemoryArray(scip, &repeat, modelData->nDays);
    SCIPallocBlockMemoryArray(scip, &numneighbors, modelData->nDays);
    for (j = 0; j < modelData->nDays; j++)
    {
        SCIPallocBlockMemoryArray(scip, &(isneighbor[j]), modelData->nC - 1);
        SCIPallocBlockMemoryArray(scip, &(values[j]), modelData->nC - 1);
    }
    for(i=0; i<modelData->nC-1; i++) /* set up neighborhood for each customer i */
    {
        for (j = 0; j < modelData->nDays; j++)
        {
            repeat[j] = FALSE;
            numneighbors[j] = 0;
            for(k = 0; k < modelData->nC - 1; k++)
            {
                isneighbor[j][k] = FALSE;
                values[j][k] = -1;
            }
        }
        tw1 = modelData->timeWindows[i];
        while(tw1 != NULL) /* we get a set of neighbors for each time window */
        {
            for(j=0; j<modelData->nC-1; j++) /* find potential neighbors induced by time window tw1 */
            {
                if(i!=j) /* no self loops */
                {
                    if (!isneighbor[tw1->day][j])
                    {
                        tw2 = modelData->timeWindows[j];
                        while(tw2 != NULL) /* iterate over the time windows of customer j */
                        {
                            if (tw1->day == tw2->day) /* potential neighbors need to be available on the say day */
                            {
                                /* calculate the earliest arrival time at customer j when starting service at
                                 * the upper bound of time window tw1 at customer i */
                                t1 = tw1->end_t + modelData->t_service[i] + modelData->t_travel[i][j];

                                /* if heuristic is activated, only close time windows can induce neighborhood */
                                if (maxWaitingTime && tw2->start_t - t1 >= maxWaitingTime)/* if the earliest arrival time is at least one hour
                                                                            * later than the start of time window tw2 we skip it */
                                {
                                    tw2 = tw2->next;
                                    continue;
                                }

                                /* calculate the earliest arrival time at customer j when starting service at
                                 * the lower bound of time window tw1 at customer i */
                                t1 = tw1->start_t + modelData->t_service[i] + modelData->t_travel[i][j];

                                if (t1 < tw2->end_t) /* the arrival time needs to be lower than the end of time window tw2 */
                                {
                                    isneighbor[tw1->day][j] = TRUE;
                                    sortvalue = tw2->start_t + modelData->t_travel[i][j] + modelData->t_travel_maxDev[i][j];
                                    nlast = NULL;
                                    n->id = j;
                                    if(numneighbors[tw1->day] == 0)
                                    {
                                        modelData->neighbors[i][tw1->day] = n;
                                        values[tw1->day][numneighbors[tw1->day]++] = sortvalue;
                                        n->next = NULL;
                                        SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                                        break;
                                    }

                                    for(k = 0; k < numneighbors[tw1->day]; k++)
                                    {
                                        if(sortvalue < values[tw1->day][k])
                                        {
                                            if(k == 0)
                                            {
                                                n->next = modelData->neighbors[i][tw1->day];
                                                modelData->neighbors[i][tw1->day] = n;
                                            }else{
                                                assert(nlast != NULL);
                                                n->next = nlast->next;
                                                nlast->next = n;
                                            }
                                            for(l = numneighbors[tw1->day]; l > k; l--)
                                            {
                                                values[tw1->day][l] = values[tw1->day][l - 1];
                                            }
                                            values[tw1->day][k] = sortvalue;
                                            numneighbors[tw1->day]++;
                                            SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                                            break;
                                        }
                                        if(k == 0){
                                            nlast = modelData->neighbors[i][tw1->day];
                                        }else{
                                            nlast = nlast->next;
                                        }
                                        if(k == numneighbors[tw1->day] - 1)
                                        {
                                            nlast->next = n;
                                            n->next = NULL;
                                            values[tw1->day][numneighbors[tw1->day]++] = sortvalue;
                                            SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                                            break;
                                        }
                                    }
                                    break;
                                }
                            }
                            tw2 = tw2->next;
                        }
                    }
                }
            }
            /* add the depot as final neighbor */
            if(!repeat[tw1->day])
            {
                n->id = modelData->nC - 1;
                n->next = NULL;
                nlast = modelData->neighbors[i][tw1->day];
                if(nlast == NULL)
                {
                    modelData->neighbors[i][tw1->day] = n;
                }else{
                    while (nlast->next != NULL)
                    {
                        nlast = nlast->next;
                    }
                    assert(nlast != NULL);
                    nlast->next = n;
                }
                SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                repeat[tw1->day] = TRUE;
            }
            tw1 = tw1->next;
        }
    }

    /* set up the neighborhood of the depot */
    i=modelData->nC-1;
    for(day=0; day<modelData->nDays; day++)
    {
        modelData->day_sizes[day] = 0;
        numneighbors[day] = 0;
        modelData->neighbors[i][day] = NULL;
        for (j = 0; j < modelData->nC - 1; j++)
        {
            /* customer j with a tw on day k is a neighbor of the depot on that particular day */
            tw2 = modelData->timeWindows[j];
            while (tw2 != NULL) /* iterate over the time windows of customer j */
            {
                if (day == tw2->day)
                {
                    modelData->day_sizes[day]++;
                    sortvalue = tw2->start_t + modelData->t_travel[i][j] + modelData->t_travel_maxDev[i][j];
                    nlast = NULL;
                    n->id = j;
                    if(numneighbors[day] == 0)
                    {
                        modelData->neighbors[i][day] = n;
                        values[day][numneighbors[day]++] = sortvalue;
                        n->next = NULL;
                        SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                        break;
                    }

                    for(k = 0; k < numneighbors[day]; k++)
                    {
                        if(sortvalue < values[day][k])
                        {
                            if(k == 0)
                            {
                                n->next = modelData->neighbors[i][day];
                                modelData->neighbors[i][day] = n;
                            }else{
                                assert(nlast != NULL);
                                n->next = nlast->next;
                                nlast->next = n;
                            }
                            for(l = numneighbors[day]; l > k; l--)
                            {
                                values[day][l] = values[day][l - 1];
                            }
                            values[day][k] = sortvalue;
                            numneighbors[day]++;
                            SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                            break;
                        }
                        if(k == 0){
                            nlast = modelData->neighbors[i][day];
                        }else{
                            nlast = nlast->next;
                        }
                        if(k == numneighbors[day] - 1)
                        {
                            nlast->next = n;
                            n->next = NULL;
                            values[day][numneighbors[day]++] = sortvalue;
                            SCIP_CALL( SCIPallocBlockMemory(scip, &n) );
                            break;
                        }
                    }
                    break;
                }
                tw2 = tw2->next;
            }
        }
        SCIPfreeBlockMemoryArray(scip, &(values[day]), modelData->nC - 1);
        SCIPfreeBlockMemoryArray(scip, &(isneighbor[day]), modelData->nC - 1);
    }
    SCIPfreeBlockMemoryArray(scip, &values, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &isneighbor, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &repeat, modelData->nDays);
    SCIPfreeBlockMemoryArray(scip, &numneighbors, modelData->nDays);
    SCIPfreeBlockMemory(scip, &n);
    return SCIP_OKAY;
}

SCIP_RETCODE customersNeighborhood(
        SCIP* scip,
        model_data* modelData,
        int customer
        ){
    int j,k;
    modelWindow* tw1 = NULL;
    modelWindow* tw2 = NULL;
    neighbor* nb = NULL;
    int** isneighbor;
    /* isneighbor[day][i] = 0 if no arc between customer and i has been saved
     *                    = 1 if only the arc (customer, i) has been saved
     *                    = 2 if only the arc (i, customer) has been saved
     *                    = 3 if both arcs (customer, i) and (i, customer) have been saved */

    /* alloc for first neighbor */
    SCIP_CALL( SCIPallocBlockMemory(scip, &nb) );
    SCIPallocBlockMemoryArray(scip, &isneighbor, modelData->nDays);
    for (j = 0; j < modelData->nDays; j++)
    {
        assert(modelData->neighbors[customer][j] == NULL);

        SCIPallocBlockMemoryArray(scip, &(isneighbor[j]), modelData->nC);
        for(k = 0; k < modelData->nC; k++)
        {
            isneighbor[j][k] = 0;
        }
    }

    tw1 = modelData->timeWindows[customer];
    while (tw1 != NULL)
    {
        /* depot is always reachable */
        if (isneighbor[tw1->day][modelData->nC - 1] % 2 == 0)
        {
            assert(modelData->neighbors[customer][tw1->day] == NULL);
            nb->id = modelData->nC - 1;
            nb->next = NULL;
            modelData->neighbors[customer][tw1->day] = nb;
            SCIP_CALL( SCIPallocBlockMemory(scip, &nb) );
            isneighbor[tw1->day][modelData->nC - 1] += 1;
        }
        /* can reach customer from depot */
        if(!isneighbor[tw1->day][modelData->nC - 1] / 2 < 1)
        {
            nb->id = customer;
            nb->next = modelData->neighbors[modelData->nC - 1][tw1->day];
            modelData->neighbors[modelData->nC - 1][tw1->day] = nb;
            SCIP_CALL( SCIPallocBlockMemory(scip, &nb) );
            isneighbor[tw1->day][modelData->nC - 1] += 2;
        }

        for(j=0; j<modelData->nC-1; j++) /* find potential neighbors induced by time window tw1 */
        {
            if (customer != j) /* no self loops */
            {
                if (isneighbor[tw1->day][j] < 3)
                {
                    tw2 = modelData->timeWindows[j];
                    while (tw2 != NULL) /* iterate over the time windows of customer j */
                    {
                        if (tw1->day == tw2->day) /* potential neighbors need to be available on the say day */
                        {
                            /* check if j is reachable from new customer on this day */
                            /* the arrival time needs to be lower than the end of time window tw2 */
                            if (isneighbor[tw1->day][j] % 2 == 0 && tw1->start_t + modelData->t_service[customer] + modelData->t_travel[customer][j] < tw2->end_t)
                            {
                                nb->id = j;
                                nb->next = modelData->neighbors[customer][tw1->day];
                                modelData->neighbors[customer][tw1->day] = nb;
                                SCIP_CALL( SCIPallocBlockMemory(scip, &nb) );
                                isneighbor[tw1->day][j] += 1;
                            }
                            /* check if new customer is reachable from j on this day */
                            if (isneighbor[tw1->day][j] / 2 < 1 && tw2->start_t + modelData->t_service[j] + modelData->t_travel[j][customer] < tw1->end_t)
                            {
                                nb->id = customer;
                                nb->next = modelData->neighbors[j][tw2->day];
                                modelData->neighbors[j][tw2->day] = nb;
                                SCIP_CALL( SCIPallocBlockMemory(scip, &nb) );
                                isneighbor[tw1->day][j] += 2;
                            }
                        }
                        tw2 = tw2->next;
                    }
                }
            }
        }
        tw1 = tw1->next;
    }
    for (j = 0; j < modelData->nDays; j++)
    {
        SCIPfreeBlockMemoryArray(scip, &(isneighbor[j]), modelData->nC);
    }
    SCIPfreeBlockMemoryArray(scip, &isneighbor, modelData->nDays);
    SCIPfreeBlockMemory(scip, &nb);
    return SCIP_OKAY;
}

/** detects customers without any valid time windows and deletes them*/
int deleteEmptyCustomers(
        SCIP* scip,
        model_data* modelData
        )
{
    int i, del;
    i = 0;
    del = 0; /* number of deleted customers */
    while(i < modelData->nC-1)
    {
        if(modelData->timeWindows[i] == NULL) /* check if no time available */
        {
            deleteCustomer(scip, i, modelData);
            i--; /* Now there is a new non-checked customer on the current position */
            del++;
        }
        i++;
    }
    return del;
}

SCIP_RETCODE reduceNumCustomers(
        SCIP* scip,
        model_data* modelData,
        int reducednC
        )
{

    replaceCustomer(scip, reducednC, modelData->nC - 1, modelData);
    modelData->nC = reducednC + 1;

    return SCIP_OKAY;
}

/** delete customer i from model data by overwriting him with the last customer and then overwrite
 * the last customer with the depot-node -> decrease nC */
SCIP_RETCODE deleteCustomer(
        SCIP* scip,
        int i,
        model_data* modelData
        )
{
    /* replace customer i with the last customer */
    replaceCustomer(scip, i, modelData->nC - 2, modelData);
    /* set depot one index lower */
    replaceCustomer(scip, modelData->nC-2, modelData->nC-1, modelData);

    modelData->nC--;

    return SCIP_OKAY;
}

/** replaces customer old by customer new */
SCIP_RETCODE replaceCustomer(
        SCIP* scip,
        int old,
        int new,
        model_data* modelData
        )
{
    int j,k;
    neighbor* nb = NULL;
    neighbor* lastnb = NULL;
    modelWindow* tmptw = NULL;
    if(old != new) /* copy replace data of customer old by data of customer old */
    {
        /* set pointer */
        modelData->customerIDs[old] = modelData->customerIDs[new];
        modelData->obj[old]=modelData->obj[new];
        modelData->t_service[old]=modelData->t_service[new];
        modelData->t_service_maxDev[old]=modelData->t_service_maxDev[new];
        /* time windows */
        tmptw = modelData->timeWindows[old];
        modelData->timeWindows[old] = modelData->timeWindows[new];
        modelData->timeWindows[new] = tmptw;
        /* travel times */
        modelData->t_travel[old]=modelData->t_travel[new];
        modelData->t_travelAM[old]=modelData->t_travelAM[new];
        modelData->t_travelNoon[old]=modelData->t_travelNoon[new];
        modelData->t_travelPM[old]=modelData->t_travelPM[new];
        modelData->t_travel_maxDev[old]=modelData->t_travel_maxDev[new];
        for(j=0;j<modelData->nC;j++) /* change travel times */
        {
            modelData->t_travel[j][old]=modelData->t_travel[j][new];
            modelData->t_travelAM[j][old]=modelData->t_travelAM[j][new];
            modelData->t_travelNoon[j][old]=modelData->t_travelNoon[j][new];
            modelData->t_travelPM[j][old]=modelData->t_travelPM[j][new];
            modelData->t_travel_maxDev[j][old]=modelData->t_travel_maxDev[j][new];
        }

        /* process neighborhoods */
        if(modelData->neighbors != NULL)
        {
            for(j=0; j<modelData->nDays; j++)
            {
                /* delete customer old from all neighbor lists */
                if(modelData->neighbors[old] != NULL)
                {
                    if(modelData->neighbors[old][j] != NULL) /* only look at active days of customer old */
                    {
                        for(k=0; k<modelData->nC; k++)
                        {
                            if (k != old)
                            {
                                lastnb = NULL;
                                if (modelData->neighbors[k][j] != NULL) /* search for customers that are active on day j */
                                {
                                    nb = modelData->neighbors[k][j];
                                    while(nb != NULL) /* check if customer old is neighbor of k on day j */
                                    {
                                        if (nb->id == old)
                                        {
                                            /* delete neighbor old from the list */
                                            if (lastnb == NULL) /* if customer old is at the start of the list */
                                            {
                                                modelData->neighbors[k][j] = nb->next;
                                                SCIPfreeBlockMemory(scip, &nb );
                                            }
                                            else
                                            {
                                                lastnb->next = nb->next;
                                                SCIPfreeBlockMemory(scip, &nb );
                                            }
                                            break;
                                        }
                                        lastnb = nb;
                                        nb = nb->next;
                                    }
                                }
                            }
                        }
                    }
                }
                /* change the neighbor-id of customer new */
                if(modelData->neighbors[new][j] != NULL) /* only look at active days of customer new */
                {
                    for(k=0; k<modelData->nC; k++) /* check for occurrences of neighbor new for every customer */
                    {
                        if (k != old)
                        {
                            if (modelData->neighbors[k][j] != NULL) /* search for customers that are active on day j */
                            {
                                nb = modelData->neighbors[k][j];
                                while(nb != NULL)
                                {
                                    if (nb->id == new) /* change id if customer new is found */
                                    {
                                        nb->id = old;
                                        break;
                                    }
                                    nb = nb->next;
                                }
                            }
                        }
                    }
                }
            }
            if(modelData->neighbors[old] != NULL)
            {
               SCIPfreeBlockMemoryArray(scip, &(modelData->neighbors[old]), modelData->nDays);
            }
            /* set pointer of customer old to customer new and set the other pointer to NULL */
            modelData->neighbors[old] = modelData->neighbors[new];
            modelData->neighbors[new] = NULL;
        }
    }

    return SCIP_OKAY;
}
