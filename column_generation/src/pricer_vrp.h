/**@file   pricer_vrp.h
 * @brief  pricer implementation with an exact labeling algorithmn
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __PRICER_VRP__
#define __PRICER_VRP__

#include "scip/scip.h"
#include "tools_data.h"
#include "vardata_vrp.h"
#include "labellist_vrp.h"
#include "label_vrp.h"

/*
 * Data structures
 */

/** @brief Variable pricer data used in the pricer_vrp.c "pricer" */
struct SCIP_PricerData
{
   SCIP_CONS**           conss;              /**< @param set covering constraints for the customers */
   int                   nconss;             /**< number of conss */
   int                   nC;                 /**< number of customers + depot */
   int                   nDays;              /**< number of days */
   model_data*           modeldata;          /**< model data */
   SCIP_CONSHDLR*        conshdlr;           /**< constraint handler for arc flow constraints */
   neighbor***           neighbors;          /**< local neighborhood */
   int                   lastnodeid;         /**< branching node id of last iteration */
   SCIP_Bool*            toDepot;            /**< is arc to depot still active during branching */
   SCIP_Bool**           isForbidden;        /**< matrix that indicates if an arc between two customers if forbidden due to arc flow branching */
   SCIP_Real             lastLPVal;          /**< the optimal value of the local LP */
};

/** creates the vrp variable pricer and includes it in SCIP */
extern
SCIP_RETCODE SCIPincludePricerVrp(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** added problem specific data to pricer and activates pricer */
extern
SCIP_RETCODE SCIPpricerVrpActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           conss,              /**< set covering constraints for the items */         
   int                   nconss,             /**< number of items to be packed */
   model_data*           modeldata           /**< modeldata */
   );

SCIP_RETCODE getDualValues(
   SCIP*                 scip,
   double*               dualvalues,
   SCIP_Bool             isFarkas
   );

/**
 * @return minimum reduced cost that could be collected */
double sumOfPossibleDualvalues(
   SCIP*                 scip,
   model_data*           modeldata,
   double*               dualvalues,
   int                   day,
   SCIP_Bool             isFarkas
   );

#endif
