/**@file   labellist_vrp.h
 * @brief  sorted list with labels, so that labels are picked in specified order (e.g. smallest reduced costs)
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __LABELLIST_VRP_H__
#define __LABELLIST_VRP_H__

#include "scip/scip.h"
#include "label_vrp.h"

typedef struct _label_list {
    labelVrp*                label;
    double                   value;
    struct _label_list*      next;
    struct _label_list*      prev;
    struct _label_list*      child;
    struct _label_list*      nextSibling;
    struct _label_list*      prevSibling;
    struct _label_list*      parent;
    SCIP_Bool                isPropagated;
} label_list;

/** Free labellist data */
SCIP_RETCODE labellistFree(
        SCIP*                   scip,
        label_list**            list
);

/** Create empty labellist*/
SCIP_RETCODE labellistCreate(
    SCIP*                   scip, 
    label_list**            list, 
    labelVrp*               label,
    double                  value
    );

/** Destroys a label list, i.e. free every attached labellist */
SCIP_RETCODE labellistDestroy(
    SCIP*                   scip,
    label_list*             list
    );

/** Inserts a label into a labellist */
SCIP_RETCODE labellistInsert(
    SCIP*                   scip,
    label_list**            list, 
    labelVrp*               label,
    double                  value
    );

/** Inserts a label into a labellist and save its address */
SCIP_RETCODE labellistInsertNew(
        SCIP*                   scip,
        label_list**            list,
        labelVrp*               label,
        label_list**            newNode,
        double                  value
);

/** Removes and gets first element of labellist */
SCIP_RETCODE labellistExtractFirst(
    SCIP*                   scip, 
    label_list**            list,
    labelVrp**              label
    );

/**
 * Checks if this labels dominates a label in the list or vice versa
 * @param scip scip instance
 * @param list list with all other labels
 * @param label new label to be checked for dominance
 * @return -1 if this label is dominated by a label in the list, 
 *          0 if there is no dominance relation to other labels,
 *          i >= 1 if this label dominates i labels already in the list */
int labellistDominanceCheck(
    SCIP*                   scip,
    label_list**            allLists,
    label_list**            list,
    labelVrp*               label,
    SCIP_Bool               usedlabels,
    int*                    deletedLabels,
    int*                    nUsedLabels
    );

void labellistPrint(
    label_list*             list
    );

void labellistPrintFirstShort(
    label_list*             list
    );    

int labellistLength(
    label_list*         list
    );

#endif