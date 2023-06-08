/**@file   cmain.c
 * @brief  Main program file for "instance generator for delay-robust VRPTW"
 * @author Andreas M. Tillmann, TU Braunschweig
 */
#include <stdio.h>
#include <string.h>
#include <scip/scip.h>
#include "tools_data.h"
#include "tools_misc.h"

/*#define SCIP_DEBUG*/ /* uncomment for more debugging output from SCIP */



/** main function: reads base-data files, creates full instance data and (optionally) writes it into a single instance-data file */
int main(
   int    argc,
   char** argv
   )
{
   SCIP* scip                  = NULL;
   char* baseDataPath          = NULL;
   char* dataSetName           = NULL;
   char* customerFile          = NULL;
   char* matrequestFile        = NULL;
   char* matresponseAMFile     = NULL;
   char* matresponseNoonFile   = NULL;
   char* matresponsePMFile     = NULL;
   char* metaFile              = NULL;
   char* outputFile            = NULL;
   instance_data* instanceData = NULL;
   model_data* modelData       = NULL;

   /* parse and check command line arguments */
   SCIP_CALL( readArguments(argc, argv, &baseDataPath, &dataSetName, &customerFile, &matrequestFile, &matresponseAMFile, &matresponseNoonFile, &matresponsePMFile, &metaFile, &outputFile) );
   assert( baseDataPath != NULL );
   assert( dataSetName != NULL );
   assert( customerFile != NULL );
   assert( matrequestFile != NULL );
   assert( matresponseAMFile != NULL );
   assert( matresponseNoonFile != NULL );
   assert( matresponsePMFile != NULL );
   assert( metaFile != NULL );

   /* initialize SCIP (currently only used here to utilize SCIP's memory management functions) */
   SCIP_CALL( SCIPcreate(&scip) );
   assert( scip != NULL );

   /* initialize instance data struct: */
   SCIP_CALL( SCIPallocBlockMemory(scip, &instanceData) );
   //SCIP_CALL( initInstanceData(scip, instanceData, -1) );
   assert( instanceData != NULL );

   /* read customer data file (discards unused data columns and extracts number of customers; which is (assumed to be) equal the
      customer numbers in matrix-request/response data)) */
   SCIP_CALL( readCustomerData(scip,customerFile,instanceData) );

   /* read matrix request data file */
   /* ... not needed at the moment: test data spans all customers listed in resp. customers.csv, so matrix-request.json
      contains no new information, because customer IDs and coordinates have already been extracted from csv-file and
      order of customers in matrix-request/response is identical to their order in csv-file ... */

   /* read matrix-response data files: morning, noon, afternoon -> three sets of distance and travel-time matrices
      (between all customers; asymmetric in general) */
   SCIP_CALL( readTravelData(scip,matresponseAMFile,matresponseNoonFile,matresponsePMFile,instanceData) );

   /* read meta data files: depot coordinates, distances and travel times (AM, noon, PM) to and from customers, start
      and end date of planning period ("default" dist./time data, similar to matrix-response.json, is currently
      ignored) */
   SCIP_CALL( readMetaData(scip,metaFile,instanceData) );



   /* create (init. and fill) model data struct: */
   SCIP_CALL( SCIPallocBlockMemory(scip, &modelData) );
   assert( modelData != NULL );
   SCIP_CALL( createModelData(scip, modelData, instanceData) );

   /* write model data struct to output file (if specified; otherwise nothing happens) */
   SCIP_CALL( writeModelData(modelData,outputFile) );


   /* de-init. model data and instance data (memory clean-up) */
   SCIP_CALL( deinitModelData(scip, modelData) );
   SCIPfreeBlockMemory(scip, &modelData);
   SCIP_CALL( deinitInstanceData(scip, instanceData) );
   SCIPfreeBlockMemory(scip, &instanceData);

   /* clean up remaining memory and exit */
   SCIP_CALL( SCIPfree(&scip) );
   BMScheckEmptyMemory();

   if( outputFile != NULL )
      free(outputFile);
   if( metaFile != NULL )
      free(metaFile);
   if( matresponsePMFile != NULL )
      free(matresponsePMFile);
   if( matresponseNoonFile != NULL )
      free(matresponseNoonFile);
   if( matresponseAMFile != NULL )
      free(matresponseAMFile);
   if( matrequestFile != NULL )
      free(matrequestFile);
   if( customerFile != NULL )
      free(customerFile);

   return 0;
}
