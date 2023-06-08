/**@file   tools_vrp.h
 * @brief  debug parameters and diverse functions for column generation
 * @author Tim Niemann, TU Braunschweig
 */

#ifndef __TOOLS_VRP__
#define __TOOLS_VRP__

#include "scip/scip.h"
#include "tools_data.h"
#include "postprocessing_vrp.h"

#define STOP_IF_FEASIBLE            FALSE       /* SCIP_BOOL,  stop the solution process when one feasible solution is found */
#define NO_REDCOST_PRICING          FALSE       /* SCIP_BOOL,  if true, there will be no reduced cost pricing, just farkas pricing to find a feasible solution */
#define PRINT_GENERATED_TOURS       FALSE       /* SCIP_BOOL,  print all generated tours to console, this does not significantly effects the runtime of SCIP*/
#define PRINT_EXACT_LABELING        FALSE       /* SCIP_BOOL,  if true, debug information for exact labeling is printed */
#define ADD_LABELS_POSITIVE_COST    FALSE       /* SCIP_BOOL,  if true, the labels with smallest, positive reduced costs are added if there are none with negative cost */
#define COMPACT_FORMULATION         FALSE       /* SCIP_BOOL,  if true, apply compact formulation pricing */

#define SAME_OBJECTIVES             TRUE        /* SCIP_BOOL,  if true, all customer dependent values for the objective function will be set to 1.0 */
#define PARALLEL_LABELING           TRUE        /* SCIP_BOOL,  if true, the exact labeling algorithm will be executed in parallel with each day as a different thread */
#define INFEASIBILITY_RECOVERY      FALSE       /* SCIP_BOOL,  if true, after detecting/estimating infeasibility, the same instance will be restarted with some customers as optional */
#define HEURISTIC_DOMINANCE         FALSE        /* SCIP_BOOL,  if true, the dominance check will be performed as a heurisitic and ignores some conditions */
#define HEURISTIC_COLLECTABLE       FALSE       /* SCIP_BOOL,  if true, the collectable reduced costs will be estimated in heuristic manner */
#define TIME_DEPENDENT_TRAVEL_TIMES TRUE        /* SCIP_BOOL,  if true, the traveltime between two customers depends on the starttime at the first customer. */

#define MIN_REQUIRED_LABELS         30          /* INT,        defines how many labels with negative reduced costs or positive farkas value must be generated before adding to master problem starts */         
#define MAX_ADDED_LABELS            1           /* INT,        defines how many labels could be added to master problem in each iteration */
#define MAX_CREATED_LABELS          1000        /* INT,        upper bound for the number of propagation steps before labeling is cancelled */
#define RELATIVE_GAP_LIMIT          0.00        /* DOUBLE,     solving stops if the relative gap is below this limit */
#define LABELING_TIME_LIMIT         30          /* INT,        time limit in seconds, after which one exact pricing iteration stops */
#define SOLVING_TIME_LIMIT          3600        /* INT,        time limit in seconds, after which scip stops */

#define MAX_HEURISTIC_TOURS         3           /* INT,        maximum number of new columns for the local search pricing */
#define WINDOW_WEIGHT_FACTOR        100         /* INT,        a fix factor for the window weight term */
#define PRICE_COLLECTING_WEIGHT     1000        /* INT,        a fix factor for the price collecting term */
#define MAX_GREEDY_TIME             600          /* INT,        maximum time that is used for the inital greedy algorithm */

#define NOON_START                  39600       /* INT,        time when AM is over and noon starts */
#define NOON_END                    54000       /* INT,        time when noon is over and PM starts */
#define EVENING_START               64800       /* INT,        time when PM is over and evening starts */
#define WORKTIME_LIMIT              39600       /* INT,        maximum worktime in seconds (31680 sec = 8h 48 min = 8h + 10%, 39600 sec = 11h) */

#define INT_BIT_SIZE                ( (int) sizeof(int) * 8 )
#define SetBit(A,k)                 ( A[(k)/INT_BIT_SIZE] |= (1 << ((k)%INT_BIT_SIZE)) )
#define ClearBit(A,k)               ( A[(k)/INT_BIT_SIZE] &= ~(1 << ((k)%INT_BIT_SIZE)) )
#define TestBit(A,k)                ( A[(k)/INT_BIT_SIZE] & (1 << ((k)%INT_BIT_SIZE)) )
#define IsSubsetBit(i,j)            ( ( ~(i) | (j) ) == -1 )


typedef struct _tuple
{
    float value;
    int index;
} tuple;

/**
 * compare function for qsort with a tuple structure
 */
int cmp_vrp(
    const void *a, 
    const void *b
    );

/**
 * Checks if the given node is the depot node 
 * @param modeldata underlying model data which contains the depot node
 * @param node node to check 
 * @return TRUE or FALSE */
SCIP_Bool nodeIsDepot(
    model_data*     modeldata,
    int             node
    );

/**
 * Computes the daytime dependent traveltime from one customer to the next one
 * @param modeldata model data with all traveltimes
 * @param start start node
 * @param end end node
 * @param starttime starttime of travelling, in most cases endtime of service, this time is used to get the actual traveltime
 * @return daytime dependent traveltime from start to end */
int getTravelTime(
    model_data*     modeldata,
    int             start,
    int             end,
    int             starttime
    );

/**
 * Find maximum between two or more integer variables
 * @param args Total number of integers
 * @param ... List of integer variables to find maximum
 * @return Maximum among all integers passed
 */
int max(
    int args, 
    ...
    );

/** Gives the next possible time window for a service start at customer 'node' on day 'day' dependent on the current time
 * @param modeldata underlying model data
 * @param node given customer
 * @param day day for the time windows
 * @param arrivaltime arrival time without delays
 * @return time window of customer 'node' with lower bound = min { tw_lb | time window tw with tw_ub >= arrival time } */
modelWindow* getNextTimeWindow(
        model_data*     modelData,
        int             node,
        int             day,
        int             arrivaltime
);

/** Gives the worst-case arrival time at node end,
 * when starting at node start and considering the current robust arrival times
 * @param modeldata underlying modeldata, contains the graph
 * @param timewindow the time window of customer end
 * @param robusttime robust arrival times of node start
 * @param start start node of the edge
 * @param end end node of the edge
 * @return if latest arrival time <= lower bound:       lower bound of time window
 *         else if latest arrival time <= upper bound:  latest arrival time
 *         else:                                        -1 */
int getLatestArrival(
        model_data*     modelData,
        modelWindow*    timewindow,
        int*            robusttime,
        int             start,
        int             end
);

/** @return Number of neighbors of node on a certain day */
int getNumberOfNeighbors(
        SCIP*           scip,
        model_data*     modelData,
        int             node,
        int             day
);

SCIP_RETCODE sortNeighborsOfNode(
        SCIP*           scip,
        model_data*     modelData,
        int             day,
        int             node,
        tuple*          list,
        int            nelements,
        double*         values
);

/** @return Number of days on which node could be served */
int getNumActiveDaysOfCustomer(
        model_data*     modelData,
        int             node
);

SCIP_RETCODE rearrangeTour(
        SCIP*                 scip,
        model_data*           modelData,
        int*                  tour,
        int                   tourlength,
        double*               tourobj,
        int                   day
);

/** Adds a node to a tour in the cheapest way.
 * @return returns FALSE if node could not be added */
SCIP_Bool addNodeToTour(
        SCIP*           scip,
        model_data*     modelData,
        int*            tour,
        int             *tourlength,
        double          *tourobj,
        int             node,
        int             day,
        int             *newpos,
        double          threshold
);

SCIP_RETCODE localSearch(
        SCIP*               scip,
        model_data*         modelData,
        int**               tours,
        int*                tourlength,
        double*             tourobj,
        int*                dayofnode,
        int                 maxdepth
);

SCIP_Bool isDeleteAllowed(
        int*                tour,
        int                 length,
        int                 depot,
        int                 pos,
        SCIP_Bool**         isForbidden
);

SCIP_Bool isExchangeAllowed(
        int*                tour,
        int                 length,
        int                 depot,
        int                 pos,
        int                 new,
SCIP_Bool**         isForbidden
);

SCIP_RETCODE twoNodeShift(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode
);

SCIP_RETCODE shiftNodes(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode
);

SCIP_RETCODE addUnvisitedNodes(
        SCIP*                 scip,
        model_data*           modelData,
        int**                 tours,
        int*                  tourlength,
        double*               tourobj,
        int*                  dayofnode,
        SCIP_Bool             isGreedy
);

/** Converts the robust arrival times of customer start
 * into the robust arrival times of customer end by using the arc (start, end) */
extern
SCIP_RETCODE updateRobustTimes(
        model_data*     modelData,
        int*            robusttime,
        int             start,
        int             end,
        int             day,
        SCIP_Bool       *isfeasible,
        double          *windowweight,
        SCIP_Bool       gethardtw,
        modelWindow**    hardwindow
);


double computeObjValue(
        SCIP*                 scip,
        model_data*           modeldata,           /**< model data */
        solutionWindow***     solutionwindows,
        SCIP_Bool*            isfeasible,
        int*                  tour,
        int*                  duration,
        int                   tourlength,
        int                   day
);

int getMSTcost(
        SCIP*                 scip,
        model_data*           modeldata
);

/**
 * Computes a lower bound for the master problem based on simple MST-costs.
 **/
SCIP_RETCODE computeStaticDualbound(
   SCIP*                 scip,               /**< SCIP data structure */
   model_data*           modeldata           /**< model data */
);

#endif
