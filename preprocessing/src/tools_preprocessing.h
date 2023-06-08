/**@file   tools_preprocessing.h
 * @brief  functions for handling of the preprocessing of the instances
 * @author Lukas Schürmann University Bonn
 */

#ifndef __TOOLS_PREPROCESSING_H__
#define __TOOLS_PREPROCESSING_H__

#include <scip/scip.h>
#include "tools_data.h"

/* read comand line arguments */
extern
SCIP_RETCODE readArguments(int argc, char** argv, char** dataFile, char** outputFile, int *maxWaitingTime, int *maxDays, int *maxTime, int *reducednC, double *windowPercentage);

/* preprocess the neighborhoods */
extern
SCIP_RETCODE processNeighbors(SCIP* scip, model_data* modelData, int deletedCustomer, int heuristic);

extern
SCIP_RETCODE customersNeighborhood(SCIP* scip, model_data* modelData, int customer);

/* merges all pairs of successive time windows */
extern
SCIP_RETCODE mergeTimeWindows(SCIP* scip, model_data* modelData);

/* delete customers with not time windows */
extern
int deleteEmptyCustomers(SCIP* scip, model_data* modelData);

extern
SCIP_RETCODE reduceNumCustomers(SCIP* scip, model_data* modelData, int reducednC);

/* deletes customer i from model data */
extern
SCIP_RETCODE deleteCustomer(SCIP* scip, int i, model_data* modelData);

/* replace customer i by customer j */
extern
SCIP_RETCODE replaceCustomer(SCIP* scip, int i, int j, model_data* modelData);

/* creates a time window from shift start to shift end for each day for each customer */
extern
SCIP_RETCODE fullFlexibility(SCIP* scip, model_data* modelData);

/** removes a given percentage of all time windows, randomly chosen. 
 * The chance of removal is applied for every window independently, which creates a binomial distribution. */
SCIP_RETCODE removeTimeWindowsBinomial(SCIP* scip, model_data* modelData, double percentage);

/** removes a given percentage of all time windows, randomly chosen.
 *  The chance of removal is applied for all windows of a customer together. */
SCIP_RETCODE removeTimeWindowsEvenly(SCIP* scip, model_data* modelData, double percentage);

/* deletes/cuts every time window that starts/ends after value maxTime */
extern
SCIP_RETCODE cutTimeWindows(SCIP* scip, model_data* modelData, int maxTime);

/* deletes time windows of days after maxDay */
extern
SCIP_RETCODE deleteDays(SCIP* scip, model_data* modelData, int maxDays);


#endif
