/**@file   cmain.c
 * @brief  Main program file for "evaluation for delay-robust VRPTW"
 * @author Lukas Schuermann, Universität Bonn
 */
#include <stdio.h>
#include <string.h>
#include <scip/scip.h>
#include "tools_evaluating.h"
#include "tools_data.h"


/*#define SCIP_DEBUG*/ /* uncomment for more debugging output from SCIP */

/** read command line arguments */
static
SCIP_RETCODE readArguments(
        int          argc,                /**< number of shell parameters */
        char**       argv,                /**< array with shell parameters */
        char**       dataFile,            /**< path of data file (e.g., ../inst_bayern.dat) */
        char**       solutionFile,        /**< path of solution file */
        char**       outputFile,          /**< optional: output instance file (incl. path; default: inst_<dataSetName>.dat) */
        double*      devfac,              /**< optional: factor for the deviation of service and travel time */
        int*         numsims,             /**< number of simulations */
        SCIP_Bool*   verbose              /**< print results in console */
)
{
    char usage[SCIP_MAXSTRLEN];
    char strbuffer[SCIP_MAXSTRLEN];
    char* locstr;
    int status;

    assert( argc >= 2 );
    assert( argv != NULL );
    assert( dataFile != NULL );
    assert( outputFile != NULL );

    /* init usage text */
    status = snprintf(usage, SCIP_MAXSTRLEN - 1,
                      "usage: %s <path of file> \
                      [-o <output instance file (optional; default used if no name is provided or if it begins with \"-\")>] \
                      [-d <deviation factor (mandatory; default = 1)>] \
                      [-n <number of simulations (mandatory; default = 100)>]", argv[0]);
    assert( 0 <= status && status < SCIP_MAXSTRLEN );

    /* init arguments */
    *dataFile           = NULL;
    *solutionFile       = NULL;
    *outputFile         = NULL;

    /* first mandatory argument: path to base data folder */
    *dataFile = argv[1];
    *solutionFile = argv[2];
    if ( *dataFile == NULL )
    {
        fprintf(stderr, "No path of model data supplied.\n");
        fprintf(stderr, "%s\n", usage);
        return SCIP_ERROR;
    }
    if ( *solutionFile == NULL )
    {
        fprintf(stderr, "No path of solution data supplied.\n");
        fprintf(stderr, "%s\n", usage);
        return SCIP_ERROR;
    }

    /* check for optional arguments */
    for (int i = 2; i < argc; i++)
    {
        /* output file */
        if ( ! strcmp(argv[i], "-o") )
        {
            assert( *outputFile == NULL );
            *outputFile = (char *) malloc ( SCIP_MAXSTRLEN * sizeof(char) +1 );

            if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
            {
                sprintf(strbuffer, "evaluation.json");
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
        /* deviation factor */
        if( ! strcmp(argv[i], "-d") )
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing deviation factor. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            sscanf(locstr, "%lf", devfac);
            if(*devfac < 0)
            {
                fprintf(stderr, "Invalid deviation factor -> Choose from interval [0, infinity). ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-n"))
        {
            if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
            {
                fprintf(stderr, "Missing number of simulations. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            i++;
            locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
            strcpy(locstr, argv[i]);
            *numsims = atoi(locstr);
            if(*numsims < 1 || *numsims > 1000)
            {
                fprintf(stderr, "Invalid number of simulations -> Choose from interval [1, 1000]. ");
                SCIPerrorMessage("%s\n", usage);
                return SCIP_ERROR;
            }
            free(locstr);
        }
        if ( ! strcmp(argv[i], "-p"))
        {
            *verbose = TRUE;
        }
    }
    return SCIP_OKAY;
}

/** main function: reads base-data files, creates full instance data and (optionally) writes it into a single instance-data file */
int main(
   int    argc,
   char** argv
   )
{
   SCIP* scip                  = NULL;
   char* dataFile              = NULL;
   char* solutionFile          = NULL;
   char* outputFile            = NULL;
   model_data* modelData       = NULL;
   solution_data* solutionData = NULL;
   double devfac               = 1.0;
   int numsims                 = 100;
   SCIP_Bool verbose           = FALSE;


   /* parse and check command line arguments */
   SCIP_CALL( readArguments(argc, argv, &dataFile, &solutionFile, &outputFile, &devfac, &numsims, &verbose) );
   assert( dataFile != NULL );

    /* initialize SCIP (currently only used here to utilize SCIP's memory management functions) */
    SCIP_CALL( SCIPcreate(&scip) );
    assert(scip != NULL);

    /* initialize model data struct: */
    SCIP_CALL( SCIPallocBlockMemory(scip, &modelData) );
    assert( modelData != NULL );

    /* initialize solution data struct */
    SCIP_CALL( SCIPallocBlockMemory(scip, &solutionData) );
    assert( solutionData != NULL );

    /* read model data file  */
    SCIP_CALL( readModelData(scip, dataFile, modelData) );

    SCIP_CALL( readSolutionData(scip, solutionFile, solutionData));

    SCIP_CALL( evaluateSolution(scip, modelData, solutionData, outputFile, devfac, numsims, verbose));

    SCIP_CALL( deinitModelData(scip, modelData) );
    SCIPfreeBlockMemory(scip, &modelData);

    SCIP_CALL( deinitSolutionData(scip, solutionData) );
    SCIPfreeBlockMemory(scip, &solutionData);
    /* clean up remaining memory and exit */
    SCIP_CALL( SCIPfree(&scip) );

   if( outputFile != NULL )
      free(outputFile);

   return 0;
}
