/**@file   cmain.c
 * @brief  Main program file for "preprocessing for delay-robust VRPTW"
 * @author Lukas Schuermann, Universität Bonn
 */
#include <stdio.h>
#include <string.h>
#include <scip/scip.h>
#include "tools_preprocessing.h"
#include "tools_data.h"


/*#define SCIP_DEBUG*/ /* uncomment for more debugging output from SCIP */


/** main function: reads base-data files, creates full instance data and (optionally) writes it into a single instance-data file */
int main(
   int    argc,
   char** argv
   )
{
   SCIP* scip                  = NULL;
   char* dataFile              = NULL;
   char* outputFile            = NULL;
   model_data* modelData       = NULL;
   int deletedCustomers        = 0;
   int maxWaitingTime          = 0;
   int maxDays                 = 365;
   int maxTime                 = DEFAULT_SHIFT_END;
   int reducednC               = 0;
   double windowPercentage     = 0.0;

   /* parse and check command line arguments */
   SCIP_CALL( readArguments(argc, argv, &dataFile, &outputFile, &maxWaitingTime, &maxDays, &maxTime, &reducednC, &windowPercentage) );
   assert( dataFile != NULL );

    /* initialize SCIP (currently only used here to utilize SCIP's memory management functions) */
    SCIP_CALL( SCIPcreate(&scip) );
    assert(scip != NULL);

    /* initialize model data struct: */
    SCIP_CALL( SCIPallocBlockMemory(scip, &modelData) );
    assert( modelData != NULL );

    /* read model data file  */
    SCIP_CALL( readModelData(scip, dataFile, modelData) );

    /* delete Days after maxDays if necessary
     * else check for time windows on invalid days */
    if( modelData->nDays > maxDays )
    {
        SCIP_CALL( deleteDays(scip, modelData, maxDays) );
        modelData->nDays = maxDays;
    }else{
        SCIP_CALL( deleteDays(scip, modelData, maxDays));
    }

    /* delete time windows after maxTime */
    if( maxTime < DEFAULT_SHIFT_END )
    {
        if(maxTime == 0)
        {
            SCIP_CALL(fullFlexibility(scip, modelData));
        }else {
            SCIP_CALL(cutTimeWindows(scip, modelData, maxTime));
        }
    }

    /* remove time windows randomly by the given percentage */
    if (windowPercentage > 0.0)
    {
        SCIP_CALL( removeTimeWindowsEvenly(scip, modelData, windowPercentage));
    }

    /* delete customers without time windows */
    deletedCustomers = deleteEmptyCustomers(scip, modelData);

    /* delete all customers except the first reducednC */
    if(reducednC > 0 && reducednC < modelData->nC)
    {
        deletedCustomers += modelData->nC - reducednC - 1;
        SCIP_CALL( reduceNumCustomers(scip, modelData, reducednC) );
    }
    
    /* merge successive time windows of each customer */
    SCIP_CALL( mergeTimeWindows(scip, modelData));

    /* preprocess the neighborhoods */
    if(modelData->neighbors == NULL || maxWaitingTime || deletedCustomers)
    {
        SCIP_CALL( processNeighbors(scip, modelData, deletedCustomers, maxWaitingTime) );
    }
    
    /* write the data in output file */
    SCIP_CALL( writeModelData(modelData, outputFile) );

    /* de-init. model data and instance data (memory clean-up) */
    modelData->nC += deletedCustomers;
    SCIP_CALL( deinitModelData(scip, modelData) );
    SCIPfreeBlockMemory(scip, &modelData);

    /* clean up remaining memory and exit */
    SCIP_CALL( SCIPfree(&scip) );

   if( outputFile != NULL )
      free(outputFile);

   return 0;
}
