/**@file   cmain.c
 * @brief  Main program file for "Delay-Resistant Robust Vehicle Routing with Heterogeneous Time Windows"
 * @author Tim Niemann, TU Braunschweig
 */

#include <stdio.h>
#include <time.h>
#include <string.h>
#include <stdlib.h>

#include "scip/scip.h"
#include "scip/scipshell.h"
#include "scip/scipdefplugins.h"

#include "tools_data.h"
#include "probdata_vrp.h"
#include "event_solution_vrp.h"
#include "pricer_vrp.h"
#include "tools_vrp.h"
#include "postprocessing_vrp.h"
#include "arcflow_branching.h"
#include "cons_arcflow.h"

/** read comand line arguments */
static
SCIP_RETCODE readArguments(
   int          argc,                /**< number of shell parameters */
   char**       argv,                /**< array with shell parameters */
   char**       inputFile,           /**< path of data file (e.g., inst_herne.dat) */
   char**       outputFile,          /**< output file (e.g. outputfile.json) */
   char**       solutionFile,        /**< path of solution file */
   int*         delayTolerance,      /**< delayTolerance, in other words: appointment window length */
   double*      weightObjective,      /**< weight for balancing the two objective functions */
   int*         gamma,                /**< robustness factor Gamma */
   double*      alphas                /**< objective function parameters */
   )
{
   int i;
   char usage[SCIP_MAXSTRLEN];
   char strbuffer[SCIP_MAXSTRLEN];
   char* locstr;
   int status;

   assert( argc >= 1 );
   assert( argv != NULL );
   assert( inputFile != NULL );

   /* init usage text */
   status = snprintf(usage, SCIP_MAXSTRLEN - 1,
      "usage: %s <path of inputfile> \
      [-w <appointment window length in sec (mandatory)>] \
      [-c <objective function weigth in [0,1] (0 = worktime, 1 = delay(default)) (mandatory)>] \
      [-o <output json file (optional; default used if no name is provided)>] \
      [-s <input solution json file (mandatory)>] \
      [-a <objective function parameters (delay, travel time, encounter probability) (default: 1, 0, 0) (mandatory) >\
      [-g <value for Gamma (mandatory)>]", argv[0]);
   assert( 0 <= status && status < SCIP_MAXSTRLEN );

   /* init arguments */
   *inputFile = NULL;
   *outputFile = NULL;
   *solutionFile = NULL;

   /* set default alphas */
   alphas[0] = 1.0;
   alphas[1] = 0.0;
   alphas[2] = 0.0;

   /* mandatory argument: inputfile */
   *inputFile = argv[1];
   if ( *inputFile == NULL )
   {
      fprintf(stderr, "No path of data supplied.\n");
      fprintf(stderr, "%s\n", usage);
      return SCIP_ERROR;
   }

   /* check for more arguments */
   for (i = 2; i < argc; i++)
   {
      /* check for (optional) appointment window length */
      if ( ! strcmp(argv[i], "-w"))
      {
         if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
         {
            fprintf(stderr, "Missing appointment window length. ");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
         strcpy(locstr, argv[i]);
         *delayTolerance = atoi(locstr);
         if(*delayTolerance < 0)
         {
            fprintf(stderr, "Invalid appointment window length -> Choose from interval [0, infinity). ");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         free(locstr);
      }
      /* check for (optional) objective function weigth */
      if ( ! strcmp(argv[i], "-c"))
      {
         if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
         {
            fprintf(stderr, "Missing objective function weigth. ");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         i++;
         locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
         strcpy(locstr, argv[i]);
         *weightObjective = atof(locstr);
         if(*weightObjective < 0 || *weightObjective > 1)
         {
            fprintf(stderr, "Invalid objective function weigth -> Choose from interval [0, 1]. ");
            SCIPerrorMessage("%s\n", usage);
            return SCIP_ERROR;
         }
         free(locstr);
      }
      /* check for (output) output file */
      if ( ! strcmp(argv[i], "-o") )
      {
         assert( *outputFile == NULL );
         *outputFile = (char *) malloc ( (int) SCIP_MAXSTRLEN * sizeof(char) );
         if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
         {
            sprintf(strbuffer, "outputfile.json");
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
       if ( ! strcmp(argv[i], "-s") )
       {
           assert( *solutionFile == NULL );
           *solutionFile = (char *) malloc ( (int) SCIP_MAXSTRLEN * sizeof(char) );
           if ( ( i == argc-1 ) || (! strncmp(argv[i+1], "-",1)) )
           {
               fprintf(stderr, "Missing solution file. ");
               SCIPerrorMessage("%s\n", usage);
               return SCIP_ERROR;
           }
           i++;
           strcpy(*solutionFile,argv[i]);
           if ( *solutionFile == NULL )
           {
               fprintf(stderr, "No path of solution supplied.\n");
               fprintf(stderr, "%s\n", usage);
               return SCIP_ERROR;
           }
           assert( i < argc );
       }
       if ( ! strcmp(argv[i], "-g") )
       {
           if( i == argc - 1 || (! strncmp(argv[i+1], "-",1)))
           {
               fprintf(stderr, "Missing robustness factor Gamma. ");
               SCIPerrorMessage("%s\n", usage);
               return SCIP_ERROR;
           }
           i++;
           locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
           strcpy(locstr, argv[i]);
           *gamma = atoi(locstr);
           if(*gamma < 0)
           {
               fprintf(stderr, "Invalid value for Gamma. Needs to be non-negative. ");
               SCIPerrorMessage("%s\n", usage);
               return SCIP_ERROR;
           }
           free(locstr);
       }
       if ( ! strcmp(argv[i], "-a") )
       {
           if( i + 3 > argc - 1 || (! strncmp(argv[i+1], "-",1)) ||
           (! strncmp(argv[i+2], "-",1)) || (! strncmp(argv[i+3], "-",1)))
           {
               printf("i %d, argc %d\n",i,argc);
               fprintf(stderr, "Missing objective parameters alpha. ");
               SCIPerrorMessage("%s\n", usage);
               return SCIP_ERROR;
           }
           for (int j=0; j < 3; j++)
           {
               i++;
               locstr = (char *) malloc ( (int)strlen(argv[i]) * sizeof(char) +1 );
               strcpy(locstr, argv[i]);
               sscanf(locstr, "%lf", &alphas[j]);
               if(alphas[j] < 0)
               {
                   fprintf(stderr, "Invalid value for alpha_%d. Needs to be non-negative. ", j);
                   SCIPerrorMessage("%s\n", usage);
                   return SCIP_ERROR;
               }
               free(locstr);
           }

       }
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE infeasibilityAnalysis(
   SCIP*                scip,
   model_data*          modelData
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_Bool* optionalCustomers;
   SCIP_Bool dayfarkasHasPositive;
   SCIP_Bool dayfarkasHasNegative;
   double* farkasvalues;
   int i, day;

   assert(scip != NULL);
   assert(modelData != NULL);
   assert(SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE);

   SCIPwarningMessage(scip, "\n*** The specified problem is infeasible. There is no solution with the given parameters for this instance. *** \n\n");
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &farkasvalues, modelData->nC - 1 + modelData->nDays) );
   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   optionalCustomers = probdata->optionalCustomers;
   probdata->optionalCost = 0;

   for (i = 0; i < modelData->nC - 1; i++)
   {
      optionalCustomers[i] = FALSE;    //overwrite decisions for the optional customers. 
      assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(probdata->conss[i]) ), "setppc", 6) );
      farkasvalues[i] = SCIPgetDualfarkasSetppc(scip, probdata->conss[i]);
   }

   dayfarkasHasNegative = FALSE;
   dayfarkasHasPositive = FALSE;
   for (day = 0; day < modelData->nDays; day++)
   {
      assert( !strncmp( SCIPconshdlrGetName( SCIPconsGetHdlr(probdata->conss[modelData->nC - 1 + day]) ), "setppc", 6) );
      farkasvalues[modelData->nC - 1 + day] = SCIPgetDualfarkasSetppc(scip, probdata->conss[modelData->nC - 1 + day]);
      if (SCIPisSumPositive(scip, farkasvalues[modelData->nC - 1 + day]))
      {
         dayfarkasHasPositive = TRUE;
      }
      if (SCIPisSumNegative(scip, farkasvalues[modelData->nC - 1 + day]))
      {
         dayfarkasHasNegative = TRUE;
      }
   }
   if (dayfarkasHasNegative && dayfarkasHasPositive)
   {
      SCIPwarningMessage(scip, "The following infeasibility analysis has a high risk of error.\n");
   }
   printf("The problem MAY be as follows:\n");

   if (dayfarkasHasPositive)
   {
      printf("There may be days for which an empty tour has to be used because not enough customers have specified time slots on these days:\n");
   }
   for (day = 0; day < modelData->nDays; day++)
   {
      if (SCIPisSumPositive(scip, farkasvalues[modelData->nC - 1 + day]))
      {
         printf(" * Day %d, customers: ", day);
         if (modelData->neighbors != NULL)
         {
            neighbor* customer = modelData->neighbors[modelData->nC - 1][day];
            while(customer != NULL)
            {
               if (SCIPisSumNegative(scip, farkasvalues[customer->id]))
               {
                  printf("%d, ", customer->id);
               }
               customer = customer->next;
            }
            printf("*\n");
         } else {
            SCIPwarningMessage(scip, "\nInfeasibilty analysis not implemented for non preprocessed data.\n");
         }
      }
   }

   if (dayfarkasHasNegative)
   {
      printf("There are some days on which certain customers would have to be planned, but this is not possible in terms of time.\n");
   }
   for(day = 0; day < modelData->nDays; day++)
   {
      if (SCIPisSumNegative(scip, farkasvalues[modelData->nC - 1 + day]))
      {
         printf(" * Day %d, customers: ", day);
         if (modelData->neighbors != NULL)
         {
            neighbor* customer = modelData->neighbors[modelData->nC - 1][day];
            while(customer != NULL)
            {
               if (SCIPisSumPositive(scip, farkasvalues[customer->id]))
               {
                  /* this customer could be a reason for infeasibilty */
                  printf("%d, ", customer->id);
                  optionalCustomers[customer->id] = TRUE;
                  probdata->optionalCost += modelData->obj[customer->id];
               }
               customer = customer->next;
            }
            printf("*\n");
         } else {
            SCIPwarningMessage(scip, "\nInfeasibilty analysis not implemented for non preprocessed data.\n");
         }
      }
   }

   if (!dayfarkasHasPositive && !dayfarkasHasNegative)
   {
      printf(" * Problematic Customers: ");
      for (i = 0; i < modelData->nC - 1; i++)
      {
         if (!SCIPisSumZero(scip, farkasvalues[i]))
            {
               /* this customer could be a reason for infeasibilty */
               printf("%d, ", i);
               optionalCustomers[i] = TRUE;
               probdata->optionalCost += modelData->obj[i];
            }
      }
      printf("\n");
   }
   SCIPfreeBlockMemoryArray(scip, &farkasvalues, modelData->nC - 1 + modelData->nDays);
   return SCIP_OKAY;
}

static
SCIP_RETCODE printVrpSolution(
   SCIP*                scip,
   model_data*          modelData
   )
{
   SCIP_SOL* sol;
   
   assert(scip != NULL);
   assert(modelData != NULL);
   assert(SCIPgetStatus(scip) != SCIP_STATUS_INFEASIBLE);
   
   sol = SCIPgetBestSol(scip);
   if (sol != NULL)
   {
      SCIP_CALL( SCIPprintSol(scip, sol, NULL, FALSE) );
   }
   return SCIP_OKAY;
}

static
SCIP_RETCODE setUpScip(
   SCIP**               scip
   )
{
   /* initialize SCIP */
   SCIP_CALL( SCIPcreate(scip) );

   /* include arc flow branching and branching data */
   SCIP_CALL( SCIPincludeBranchruleArcFlow(*scip) );
   SCIP_CALL( SCIPincludeConshdlrArcFlow(*scip) );

   /* include vrp pricer  */
   SCIP_CALL( SCIPincludePricerVrp(*scip) );

   /* include event handler for new solutions */
   SCIP_CALL( SCIPincludeEventHdlrBestsol(*scip) );

   /* include default SCIP plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(*scip) );

   /* change display columns */
   SCIP_CALL( SCIPsetIntParam(*scip,"display/nfrac/active",2) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/maxdepth/active",0) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/curcols/active",2) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/cuts/active",2) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/lpobj/active",2) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/primalgap/active",0) );
   SCIP_CALL( SCIPsetIntParam(*scip,"display/gap/active",2) ); 
   SCIP_CALL( SCIPsetIntParam(*scip,"display/lpavgiterations/active",0) );
   SCIP_CALL( SCIPsetIntParam(*scip,"randomization/randomseedshift",1) );

    /* for column generation instances, disable restarts */
   SCIP_CALL( SCIPsetIntParam(*scip,"presolving/maxrestarts",0) );

   /* turn off all separation algorithms */
   SCIP_CALL( SCIPsetSeparating(*scip, SCIP_PARAMSETTING_OFF, TRUE) );

   /* Stop after one feasible solution, if set */
   if (STOP_IF_FEASIBLE)
   {
      SCIP_CALL( SCIPsetLongintParam( *scip, "constraints/countsols/sollimit", 1) );
   }

   /* set the limit for the relative gap after which the solving process is stopped */
   SCIP_CALL (SCIPsetRealParam(*scip, "limits/gap", RELATIVE_GAP_LIMIT));

   /* set a time limit */
   SCIP_CALL( SCIPsetRealParam(*scip, "limits/time", SOLVING_TIME_LIMIT));
   return SCIP_OKAY;
}

/** creates a SCIP instance with default plugins, evaluates command line parameters, runs SCIP appropriately,
 *  and frees the SCIP instance
 */
static
SCIP_RETCODE runColumnGenerationModel(
   int                   argc,               /**< number of shell parameters */
   char**                argv                /**< array with shell parameters */
   )
{
   SCIP* scip = NULL;
   SCIP_PROBDATA* probdata = NULL;
   model_data* modelData = NULL;
   solution_data* solutionData = NULL;
   char* inputfile = NULL;
   char* outputfile = NULL;
   char* solutionfile = NULL;
   int delayTolerance = 0;
   double weightObjective = 1.0;
   int gamma = -1;
   double* alphas;
   SCIP_Bool* optionalCustomers = NULL;
   int i;

   alphas = malloc(3*sizeof(double));

   /*********
    * Setup *
    *********/

   SCIP_CALL( setUpScip(&scip));
   assert(scip != NULL);

   SCIP_CALL( readArguments(argc, argv, &inputfile, &outputfile, &solutionfile, &delayTolerance, &weightObjective, &gamma, alphas) );
   assert(inputfile != NULL);

   /*********************
    * Read Data
    *********************/

   SCIP_CALL( SCIPallocBlockMemory(scip, &modelData) );
   assert(modelData != NULL);
   SCIP_CALL( readModelData(scip, inputfile, modelData) );
   if (modelData->neighbors == NULL)
   {
      SCIPwarningMessage(scip, "The given instance is not preprocessed! The computation could be very slow and inefficient!\n");
   }
   if(solutionfile != NULL)
   {
       SCIP_CALL( SCIPallocBlockMemory(scip, &solutionData) );
       assert( solutionData != NULL );
       printf("Start reading solutiondata");
       SCIP_CALL( readSolutionData(scip, solutionfile, solutionData) );
       printf("Finished reading solutiondata");
       free(solutionfile);
       if(solutionData->nDays > modelData->nDays || solutionData->nC > modelData->nC)
       {
           fprintf(stderr, "Solution file does not fit to model data - number of customers\n");
           return SCIP_ERROR;
       }
       if(solutionData->windowLength != delayTolerance)
       {
           fprintf(stderr, "Solution file does not fit to model data - delay tolerance (%d != %d)\n", delayTolerance, solutionData->windowLength);
           return SCIP_ERROR;
       }
       for(i = 0; i < 3; i++)
       {
           if(solutionData->alphas[i] != alphas[i])
           {
               fprintf(stderr, "Solution file does not fit to model data - alpha_%d (%f != %f)\n", i, alphas[i], solutionData->alphas[i]);
               return SCIP_ERROR;
           }
       }
       if((gamma >= 0 && solutionData->gamma != gamma) || (gamma < 0 && solutionData->gamma != 10))
       {
           fprintf(stderr, "Solution file does not fit to model data - robustness degree (%d != %d)\n", gamma, solutionData->gamma);
           return SCIP_ERROR;
       }
   }
   if(gamma >= 0)
   {
       modelData->maxDelayEvents = gamma;
   }
   /*********************
    * Create Master Problem
    *********************/
   
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &optionalCustomers, modelData->nC - 1) );
   for (i = 0; i < modelData->nC - 1; i++)
   {
      optionalCustomers[i] = FALSE;
   }
   SCIP_CALL( SCIPprobdataCreate(scip, modelData, solutionData, delayTolerance, alphas, optionalCustomers) );
   SCIPfreeBlockMemoryArray(scip, &optionalCustomers, modelData->nC - 1);

   /********************
    * Solve Problem
    ***********************/

   SCIP_CALL( SCIPsolve(scip) );
   probdata = SCIPgetProbData(scip);
   
   /********************
    * Print Solution
    ***********************/

//   SCIP_CALL( SCIPprintStatistics(scip, NULL) );

   if (SCIPgetStatus(scip) != SCIP_STATUS_INFEASIBLE) 
   {
      printVrpSolution(scip, modelData);
   }
   
   if (outputfile != NULL && SCIPgetNSols(scip) > 0)
   {
      createOutputFile(scip, modelData, outputfile);
      free(outputfile);
   }

   /********************
    * Recovery         *
    ********************/

   if (SCIPgetStatus(scip) == SCIP_STATUS_INFEASIBLE)
   {      
      /* detect, which customers led to the infeasibility */
      SCIP_CALL( infeasibilityAnalysis(scip, modelData));
   }

   /* if there is no feasible solution to the original problem, neither through heuristics or exact methods */
   if (INFEASIBILITY_RECOVERY && probdata->optionalCost > 0)
   {
      SCIP* scipRecovery = NULL;
      printf("\n * No feasible solution found. Unserved customers: ");
      for (i = 0; i < modelData->nC-1; i++)
      {
         if (probdata->optionalCustomers[i] == TRUE)
         {
            printf("%d, ", i);
         }
      }
      printf("*\n * Restarting Scip with the problematic customers as optional. *\n\n");

      SCIP_CALL( setUpScip(&scipRecovery));
      assert(scipRecovery != NULL);

      SCIP_CALL( SCIPprobdataCreate(scipRecovery, modelData, solutionData, delayTolerance, alphas, probdata->optionalCustomers) );
      SCIP_CALL( SCIPsolve(scipRecovery) );

      if (SCIPgetStatus(scipRecovery) == SCIP_STATUS_INFEASIBLE)
      {
         SCIPwarningMessage(scipRecovery, "\nInfeasibilty recovery did not lead to a feasibile solution.\n");
      } else {
         printVrpSolution(scipRecovery, modelData);
         if (outputfile != NULL && SCIPgetNSols(scipRecovery) > 0)
         {
            createOutputFile(scipRecovery, modelData, outputfile);
            free(outputfile);
         }
      }
      SCIP_CALL( SCIPfree(&scipRecovery) );
   }

   /********************
    * Deinitialization *
    ********************/
    free(alphas);

   if(solutionData != NULL)
   {
       SCIP_CALL(deinitSolutionData(scip, solutionData));
   }
   SCIP_CALL( deinitModelData(scip, modelData) );
   SCIPfreeBlockMemory(scip, &modelData);
   SCIP_CALL( SCIPfree(&scip) );

   return SCIP_OKAY;
}

int
main(
   int                        argc,
   char**                     argv
   )
{
   SCIP_RETCODE retcode;

   retcode = runColumnGenerationModel(argc, argv);
   if( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }

   return 0;
}

