/**@file   tools_misc.c
 * @brief  miscellaneous auxiliary functions
 * @author Andreas M. Tillmann, TU Braunschweig
 */

#include "tools_misc.h"
#include "tools_data.h"
#include <string.h>
#include <stdlib.h>
#include <scip/scip.h>



/** read comand line arguments */
SCIP_RETCODE readArguments(
   int          argc,                /**< number of shell parameters */
   char**       argv,                /**< array with shell parameters */
   char**       baseDataPath,        /**< path to base-data folder (e.g., ../../../test-data/base-data) */
   char**       dataSetName,         /**< name of data set (subfolder in base-data folder, e.g., herne) */
   char**       customerFile,        /**< optional: customer (csv-)data file name (default: customers.csv in <baseDataPath>/<dataSetName>/) */
   char**       matrequestFile,      /**< optional: customer IDs and coordinates (json-)data file (incl. path; default: matrix-request.json in <baseDataPath>/<dataSetName>/) */
   char**       matresponseAMFile,   /**< optional: morning distance & travel time (customer matrix) data file (incl. path; default: <baseDataPath>/<dataSetName>/matrix-response-morning.json) */
   char**       matresponseNoonFile, /**< optional: noon distance & travel time (customer matrix) data file (incl. path; default: <baseDataPath>/<dataSetName>/matrix-response-noon.json) */
   char**       matresponsePMFile,   /**< optional: afternoon distance & travel time (customer matrix) data file (incl. path; default: <baseDataPath>/<dataSetName>/matrix-response-afternoon.json) */
   char**       metaFile,            /**< optional: meta data file (incl. path; default: <baseDataPath>/<dataSetName>/meta.json) */
   char**       outputFile           /**< optional: output instance file (incl. path; default: inst_<dataSetName>.dat) */
   )
{
   int i;
   char usage[SCIP_MAXSTRLEN];
   char strbuffer[SCIP_MAXSTRLEN];
   int status;

   assert( argc >= 2 );
   assert( argv != NULL );
   assert( baseDataPath != NULL );
   assert( dataSetName != NULL );
   assert( customerFile != NULL );
   assert( matrequestFile != NULL );
   assert( matresponseAMFile != NULL );
   assert( matresponseNoonFile != NULL );
   assert( matresponsePMFile != NULL );
   assert( metaFile != NULL );
   assert( outputFile != NULL );

   /* init usage text */
   status = snprintf(usage, SCIP_MAXSTRLEN - 1,
      "usage: %s <path to base data> <data set (subfolder) name>  [-o <output instance file (optional; default used if no name is provided or if it begins with \"-\")>] [-c <customer data file>] [-mr <matrix-request file>] [-mm <matrix-response (morning) file>] [-mn <matrix-response (noon) file>] [-ma <matrix-response (afternoon) file>] [-me <meta file>]", argv[0]);
   assert( 0 <= status && status < SCIP_MAXSTRLEN );

   /* init arguments */
   *baseDataPath        = NULL;
   *dataSetName         = NULL;
   *customerFile        = NULL;
   *matrequestFile      = NULL;
   *matresponseAMFile   = NULL;
   *matresponseNoonFile = NULL;
   *matresponsePMFile   = NULL;
   *metaFile            = NULL;
   *outputFile          = NULL;

   /* first mandatory argument: path to base data folder */
   *baseDataPath = argv[1];
   if ( *baseDataPath == NULL )
   {
      fprintf(stderr, "No path to base-data supplied.\n");
      fprintf(stderr, "%s\n", usage);
      return SCIP_ERROR;
   }
   /* second mandatory argument: name of the data set (subfolder) */
   *dataSetName  = argv[2];
   if ( *dataSetName == NULL )
   {
      fprintf(stderr, "No data set (subfolder) name supplied.\n");
      fprintf(stderr, "%s\n", usage);
      return SCIP_ERROR;
   }

   /* check for optional arguments */
   for (i = 3; i < argc; i++)
   {
      /* check for customer data file */
      if ( ! strcmp(argv[i], "-c") )
      {
         assert( *customerFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing customer file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *customerFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 ); /* +1 b/c strcpy adds a \0 character */
         assert( *customerFile != NULL );
         strcpy(*customerFile,argv[i]);
         assert( i < argc );
      }

      /* check for customer ID/coordinates matrix data file */
      if ( ! strcmp(argv[i], "-mr") )
      {
         assert( *matrequestFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing matrix-request file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *matrequestFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         assert( *matrequestFile != NULL );
         strcpy(*matrequestFile,argv[i]);
         assert( i < argc );
      }

      /* check for distance/travel times (morning) matrix data file */
      if ( ! strcmp(argv[i], "-mm") )
      {
         assert( *matresponseAMFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing matrix-response-morning file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *matresponseAMFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         assert( *matresponseAMFile != NULL );
         strcpy(*matresponseAMFile,argv[i]);
         assert( i < argc );
      }

      /* check for distance/travel times (noon) matrix data file */
      if ( ! strcmp(argv[i], "-mn") )
      {
         assert( *matresponseNoonFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing matrix-response-noon file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *matresponseNoonFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         assert( *matresponseNoonFile != NULL );
         strcpy(*matresponseNoonFile,argv[i]);
         assert( i < argc );
      }

      /* check for distance/travel times (afternoon) matrix data file */
      if ( ! strcmp(argv[i], "-ma") )
      {
         assert( *matresponsePMFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing matrix-response-afternoon file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *matresponsePMFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         assert( *matresponsePMFile != NULL );
         strcpy(*matresponsePMFile,argv[i]);
         assert( i < argc );
      }

      /* check for meta data file */
      if ( ! strcmp(argv[i], "-me") )
      {
         assert( *metaFile == NULL );
         if ( i == argc-1 )
         {
            fprintf(stderr, "Missing meta data file name.\n");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         *metaFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         assert( *metaFile != NULL );
         strcpy(*metaFile,argv[i]);
         assert( i < argc );
      }

      /* check for (output) instance file */
      if ( ! strcmp(argv[i], "-o") )
      {
         assert( *outputFile == NULL );
         *outputFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
         if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
         {
            sprintf(strbuffer, "inst_%s.dat", *dataSetName);
            strcpy(*outputFile,strbuffer);
         }
         else
         {
            i++;
            strcpy(*outputFile,argv[i]);
            assert( i < argc );
         }
         assert( *outputFile != NULL );
      }

   }
   /* set default values for arguments still unspecified */
   if( *customerFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/customers.csv", *baseDataPath, *dataSetName);
      *customerFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 ); /* +1 b/c strcpy adds a \0 character */
      strcpy(*customerFile,strbuffer);
   }
   if( *matrequestFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/matrix-request.json", *baseDataPath, *dataSetName);
      *matrequestFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
      strcpy(*matrequestFile,strbuffer);
   }
   if( *matresponseAMFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/matrix-response-morning.json", *baseDataPath, *dataSetName);
      *matresponseAMFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
      strcpy(*matresponseAMFile,strbuffer);
   }
   if( *matresponseNoonFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/matrix-response-noon.json", *baseDataPath, *dataSetName);
      *matresponseNoonFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
      strcpy(*matresponseNoonFile,strbuffer);
   }
   if( *matresponsePMFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/matrix-response-afternoon.json", *baseDataPath, *dataSetName);
      *matresponsePMFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
      strcpy(*matresponsePMFile,strbuffer);
   }
   if( *metaFile == NULL )
   {
      sprintf(strbuffer, "%s/%s/meta.json", *baseDataPath, *dataSetName);
      *metaFile = (char *) malloc ( (int)sizeof(strbuffer) * sizeof(char) +1 );
      strcpy(*metaFile,strbuffer);
   }
   /* *outputFile remains NULL if option "-o" was not used; default is set earlier in case "-o" is not followed by filename */

   return SCIP_OKAY;
}





