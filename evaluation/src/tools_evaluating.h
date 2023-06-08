/**@file   tools_evaluating.c
 * @brief  functions for evaluation solutions
 * @author Lukas Schürmann University Bonn
 */

#ifndef __TOOLS_PREPROCESSING_H__
#define __TOOLS_PREPROCESSING_H__

#include <scip/scip.h>
#include "tools_data.h"

#define NOON_START                  39600       /* INT,         time when AM is over and noon starts */
#define NOON_END                    54000       /* INT,         time when noon is over and PM starts */
#define EVENING_START               64800       /* INT,         time when PM is over and evening starts */
#define WORKTIME_LIMIT              60000       /* INT,         maximum worktime in seconds (31680 sec = 8h 48 min = 8h + 10%) */

#define RAND_SERVICE                TRUE        /* SCIP_BOOL,   if true, include random deviation for the service times */
#define RAND_TRAVEL                 TRUE        /* SCIP_BOOL,   if true, include random deviation for the travel times */
#define RAND_ENCOUNTER              FALSE       /* SCIP_BOOL,   if true, we encounter the customers with a certain probability */

#define SERVICE_LB                  30          /* INT,         service time for customers that are not available */

#define RNDSEED                        1

extern
SCIP_RETCODE evaluateSolution(SCIP* scip, model_data* modelData, solution_data* solutionData, char* outputFile,
                              double devfac, int numsims, SCIP_Bool verbose);


#endif
