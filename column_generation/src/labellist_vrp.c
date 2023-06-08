/**@file   labellist_vrp.c
 * @brief  sorted list with labels, so that labels are picked in specified order (e.g. smallest reduced costs)
 * @author Tim Niemann, TU Braunschweig
 */

#include "scip/scip.h"
#include "labellist_vrp.h"
#include "label_vrp.h"

/** Free labellist data */
SCIP_RETCODE labellistFree(
    SCIP*                   scip, 
    label_list**            list
    )
{
    assert(scip != NULL);
    assert(list != NULL);
    if (*list != NULL)
    {
        SCIPfreeMemory(scip, list);
        *list = NULL;
    }

    return SCIP_OKAY;
}

/** Create empty labellist*/
SCIP_RETCODE labellistCreate(
    SCIP*                   scip, 
    label_list**            list, 
    labelVrp*               label,
    double                  value
    ) 
{
    assert(scip != NULL);
    assert(list != NULL);
    assert(*list == NULL);
    assert(label != NULL);

    SCIP_CALL( SCIPallocMemory(scip, list) );

    (*list)->label = label;
    (*list)->value = value;
    (*list)->next = NULL;
    (*list)->prev = NULL;
    (*list)->nextSibling = NULL;
    (*list)->prevSibling = NULL;
    (*list)->parent = NULL;
    (*list)->child = NULL;
    (*list)->isPropagated = FALSE;
    
    return SCIP_OKAY;
}

/** Destroys a label list, i.e. free every attached labellist */
SCIP_RETCODE labellistDestroy(
    SCIP*                   scip,
    label_list*             list
    )
{
    label_list* tmp;

    assert(scip != NULL);
    
    while(list != NULL) 
    {
        tmp = list->next;
        labelVrpFree(scip, &(list->label));
        labellistFree(scip, &list);
        list = tmp;
    }

    return SCIP_OKAY;
}

/** Inserts a label into a labellist */
SCIP_RETCODE labellistInsert(
    SCIP*                   scip,
    label_list**            list, 
    labelVrp*               label,
    double                  value
    ) 
{
    label_list* newNode;
    label_list* tmp;
    assert(label != NULL);

    if (*list == NULL)
    {
        SCIP_CALL( labellistCreate(scip, list, label, value));
        assert(*list != NULL);
        return SCIP_OKAY;
    }
    
    newNode = NULL;
    assert(*list != NULL);
    assert((*list)->label != label);
    SCIP_CALL( labellistCreate(scip, &newNode, label, value));
    assert(newNode != NULL);

    /* Insert at the start, set new head node */
    if ((*list)->value >= value)
    {
        newNode->next = *list;
        *list = newNode;
        return SCIP_OKAY;
    }
    tmp = *list;
    while (tmp->next != NULL)
    {
        if (tmp->next->value < value)
        {
            tmp = tmp->next;
        /* Insertion in the middle */
        } else {
            assert(tmp->value <= value);
            assert(value <= tmp->next->value);
            newNode->next = tmp->next;
            tmp->next = newNode;
            return SCIP_OKAY;
        }
    }
    /* Insertion at the end */
    assert(tmp->next == NULL);
    assert(tmp->value < value);

    tmp->next = newNode;
    assert((*list)->next != (*list));
    
    return SCIP_OKAY;
}

/** Inserts a label into a labellist and save its address */
SCIP_RETCODE labellistInsertNew(
        SCIP*                   scip,
        label_list**            list,
        labelVrp*               label,
        label_list**            newNode,
        double                  value
)
{
    label_list* tmp;
    assert(label != NULL);

    /* heuristic ordering, results in modified depth search with less stagnation */
    value = value + 0.1 * label->collactableRedcost;
    if (*list == NULL)
    {
        if(*newNode == NULL)
        {
            SCIP_CALL(labellistCreate(scip, list, label, value));
            *newNode = *list;
        }else{
            *list = *newNode;
        }
        assert(*list != NULL);
        return SCIP_OKAY;
    }

    assert(*list != NULL);
    assert((*list)->label != label);
    if(*newNode == NULL)
    {
        SCIP_CALL(labellistCreate(scip, newNode, label, value));
    }
    assert(*newNode != NULL);

    /* Insert at the start, set new head node */
    if ((*list)->value >= value)
    {
        (*newNode)->next = *list;
        (*list)->prev = *newNode;
        *list = *newNode;
        return SCIP_OKAY;
    }
    tmp = *list;
    while (tmp->next != NULL)
    {
        if (tmp->next->value < value)
        {
            tmp = tmp->next;
        } else {
            /* Insertion in the middle */
            assert(tmp->value <= value);
            assert(value <= tmp->next->value);
            (*newNode)->next = tmp->next;
            (*newNode)->prev = tmp;
            tmp->next->prev = *newNode;
            tmp->next = *newNode;
            return SCIP_OKAY;
        }
    }
    /* Insertion at the end */
    assert(tmp->next == NULL);
    assert(tmp->value < value);

    (*newNode)->prev = tmp;
    tmp->next = *newNode;
    assert((*list)->next != (*list));
    return SCIP_OKAY;
}

/** Removes and gets first element of labellist */
SCIP_RETCODE labellistExtractFirst(
    SCIP*                   scip, 
    label_list**            list,
    labelVrp**              label
    ) 
{
    label_list* tmp = *list;
    assert(*list != NULL);
    assert(*label == NULL);

    *label = (*list)->label;
    *list = (*list)->next;
    labellistFree(scip, &tmp);

    return SCIP_OKAY;
}

/** deletes labellist from global labellist-tree */
static
SCIP_RETCODE deleteFromFamilyTree(
        label_list*         list
        ){
    if(list->prevSibling == NULL)
    {
        if(list->nextSibling != NULL)
        {
            list->nextSibling->prevSibling = NULL;
        }
        list->parent->child = list->nextSibling;
    }else{
        list->prevSibling->nextSibling = list->nextSibling;
        if(list->nextSibling != NULL)
        {
            list->nextSibling->prevSibling = list->prevSibling;
        }
    }
    return SCIP_OKAY;
}


/** deletes labellist from customer labellist */
static
SCIP_RETCODE deleteFromLabelList(
        label_list**        allLists,
        label_list*         list
        ){
    if(list->prev == NULL)
    {
        if(list->next != NULL) {
            list->next->prev = NULL;
        }
        *allLists = list->next;
    }else{
        list->prev->next = list->next;
        if(list->next != NULL)
        {
            list->next->prev = list->prev;
        }
    }
    return SCIP_OKAY;
}

/** deletes labellist references */
static
SCIP_RETCODE deleteList(
        SCIP*               scip,
        label_list**        allLists,
        label_list*         list
        ){
    assert(list != NULL);

    deleteFromFamilyTree(list);
    deleteFromLabelList(allLists, list);

    labelVrpFree(scip, &(list->label));
    labellistFree(scip, &list);

    return SCIP_OKAY;
}

/** recursively deletes all labels that got propagated from the current one */
static
SCIP_RETCODE deleteChildren(
        SCIP*              scip,
        label_list**       activeLabels,
        label_list**       oldLabels,
        label_list*        list,
        int*               deletedLabels,
        int*               nUsedLabels
        ){
    if(list->isPropagated)
    {
        label_list* child = list->child;
        while(child != NULL)
        {
            deleteChildren(scip, activeLabels, oldLabels, child, deletedLabels, nUsedLabels);

            child = list->child;
        }
        nUsedLabels[list->label->node]--;
        deleteList(scip, &oldLabels[list->label->node], list);

        return SCIP_OKAY;
    }else{
        deleteList(scip, &activeLabels[list->label->node], list);
        (*deletedLabels)++;

        return SCIP_OKAY;
    }
}

/**
 * Checks if this labels dominates a label in the list or vice versa
 * @param scip scip instance
 * @param list list with all other labels
 * @param label new label to be checked for dominance
 * @return -1 if this label is dominated by a label in the list, 
 *          0 if there is no dominance relation to other labels,
 *          1 if this label dominates a label already in the list */
int labellistDominanceCheck(
    SCIP*                   scip,
    label_list**            activeLabels,
    label_list**            oldLabels,
    labelVrp*               label,
    SCIP_Bool               usedlabels,
    int*                    deletedLabels,
    int*                    nUsedLabels
    )
{
    label_list *tmplist;
    label_list *nextlist;
    int i = 0;

    assert(label != NULL);

    if(!usedlabels)
    {
        tmplist = activeLabels[label->node];
    }else{
        tmplist = oldLabels[label->node];
    }

    while (tmplist != NULL)
    {
        nextlist = tmplist->next;
        /* a member of the list dominates the new label */
        if (labelVrpDominates(scip, tmplist->label, label))
        {
            assert(i == 0);
            return -1;
        }
        /* or the new label dominates a member of the list, delete the already contained one */
        else if (labelVrpDominates(scip, label, tmplist->label))
        {
            if(!usedlabels)
            {
                i++;
                deleteList(scip, &activeLabels[label->node], tmplist);
            }else{
                /* delete all descendants of tmplist */
                SCIP_CALL( deleteChildren(scip, activeLabels, oldLabels, tmplist, deletedLabels, nUsedLabels) );
            }
        }
        tmplist = nextlist;
    }

    return i;
}

void labellistPrint(
    label_list*             list
    ) 
{    
    while (list != NULL)
    {
        labelVrpPrint(list->label);
        assert(list->next != list);
        list = list->next;
    }
}

void labellistPrintFirstShort(
    label_list*             list
    )
{
    int i;
    if (list != NULL)
    {
        for (i = 0; i < list->label->nvisitednodes; i++)
        {
            printf("%d->", list->label->visitednodes[i]);
        }
        printf("\n");
    }
}

int labellistLength(
    label_list*         list
    )
{
    int i = 0;
    while(list != NULL)
    {
        list = list->next;
        i++;
    }
    return i;
}