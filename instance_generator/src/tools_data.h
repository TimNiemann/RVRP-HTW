/**@file   tools_data.h
 * @brief  functions for data handling (reading, writing, processing)
 * @author Andreas M. Tillmann, TU Braunschweig
 */

#ifndef __TOOLS_DATA_H__
#define __TOOLS_DATA_H__

#include <scip/scip.h>

#define REFDATE_DAY            1
#define REFDATE_MONTH          1
#define REFDATE_YEAR           2000
#define REFDATE_WEEKDAY        5     /* <-> Saturday; we start counting weekdays at Monday <-> 0 */

#define DEFAULT_SHIFT_START    21600 /* <-> 06:00; time in seconds when driver's daily work shift starts (hard lower bound for departure from depot) */
#define DEFAULT_SHIFT_END      79200 /* <-> 22:00; time in seconds when driver's daily work shift ends   (hard upper bound for return to depot) */
#define DEFAULT_BREAK_START    39600 /* <-> 11:00; time in seconds when the driver can start its break */
#define DEFAULT_BREAK_END      54000 /* <-> 14:00; time in seconds until when the driver has to take its break */
#define DEFAULT_BREAK_LENGTH   1800  /* <-> 30 minutes; length of the break */
#define DEFAULT_MAXDELAYEVENTS 10    /* <-> maximal number of "delay events" the optimization models consider in robustness part */

/** date structure */
typedef struct _date
{
   int                   day;                    /**< day (1-31) */
   int                   month;                  /**< month (1-12) */
   int                   year;                   /**< year */
   int                   weekday;                /**< weekday index (0-6) */
}date;

/** possible delivery time windows (linked list) structure [for instance data] */
typedef struct _timeWindow
{
   date*                 twDate;                 /**< date struct; (day,month,year can be -1 if only weekday is relevant) */
   int                   start_hour;             /**< hour (24h format) when time window opens */
   int                   start_min;              /**< minute when time window opens */
   int                   end_hour;               /**< hour (24h format) when time window closes */
   int                   end_min;                /**< minute when time window closes */
   SCIP_Bool             recurring;              /**< Boolean indicator whether timeWindow is linked to specific date (FALSE) or recurring weekday (TRUE) */
   double                weigth;                 /**< weight of a time window (probability of availability of the customer at that time window) */
   struct _timeWindow*   next;                   /**< pointer to next timeWindow in list (or NULL) */
}timeWindow;

typedef struct _dayProperty
{
   int                   index;
   int                   shift_start;
   int                   shift_end;
   int                   break_start;
   int                   break_end;
   int                   break_length;
}dayProperty;

/** possible delivery time windows (linked list) structure [for model data] */
typedef struct _modelWindow
{
   int                   day;                    /**< day (index w.r.t. planning period) */
   int                   start_t;                /**< second when time window opens (relative to beginning of resp. work day) */
   int                   end_t;                  /**< second when time window closes (relative to beginning of resp. work day) */
   double                weigth;                 /**< weight of a time window (probability of availability of the customer at that time window) */
   struct _modelWindow*  next;                   /**< pointer to next timeWindow in list (or NULL) */
}modelWindow;

/** list of neighbor nodes for each customer node and day of planning period */
typedef struct _neighbor
{
   int                   id;                     /**< index of a customer */
   struct _neighbor*     next;                   /**< pointer to next neighbor in list (or NULL) */
}neighbor;

typedef struct _solutionWindow
{
    int                 day;                    /**< day (index w.r.t. planning period) */
    int                 start_t;                /**< second when time window opens (relative to beginning of resp. work day), equals earliest arrival */
    int                 end_t;                  /**< second when time window closes (relative to beginning of resp. work day), eqials earliest arrival plus given time window length */
    int                 latest_arrival;         /**< second when customer is served in the computed worst case szenario */
    modelWindow*        modelwindow;
} solutionWindow;

/** solution data (stores the computed solution, in a form that can be written to an output file) */
typedef struct _solutionData
{
    int                 nC;                     /**< number of customers + 1 (depot stored as last "customer") */
    int                 nDays;                  /**< number of work days in planning period (from startDate to endDate, including both) */
    int                 gamma;                  /**< robustness degree Gamma */
    int                 windowLength;           /**< length of the appointment windows */
    SCIP_Real           objCoefficient;         /**< parameter c in the objective function */
    solutionWindow**    timeWindows;            /**< nC-array to store computed delivery time windows for each customer */
    SCIP_Real*          obj;                    /**< nDays-array with objective function coefficients for each tour */
    int**               tour;                   /**< (nDays x tourLength)-jagged-array with the sequence of visited customers for each tour */
    int*                tourLength;             /**< nDays-array with number of visited customers */
    int*                tourDuration;           /**< nDays-array with worst case tour duration in seconds for each tour */
    int*                customerIDs;            /**< nC-array with customer IDs for output/information purposes; depot gets ID -1 */
    SCIP_Real*          alphas;
} solution_data;

/** instance data (stores base-data) */
typedef struct _instanceData
{
   int                   nC;                     /**< number of customers */
   int                   nDays;                  /**< number of work days */
   int*                  customerIDs;            /**< nC-array with customer IDs (integers) */
   SCIP_Real**           coordinates;            /**< (nC x 2)-array with customer coordinates (doubles) */
   timeWindow**          timeWindows;            /**< nC-array with linked list structures to store customers' delivery time windows */
   int*                  deliveryFreq;           /**< nC-array with customer delivery frequencies (in weeks; nonneg. integers) */
   int*                  availability;           /**< nC-array with customer availability values (binary integers) */
   SCIP_Real*            estim_availability;     /**< nC-array with customer availability probability/estimate values (doubles in [0,1]) */
   int*                  appointment_sens;       /**< nC-array with customer appointment sensibility values (1,2,3 representing characters A,B,C) */
   int*                  appointment_tol;        /**< nC-array with customer appointment (rescheduling) tolerance values (nonneg. integers in [0,20]) */
   int*                  t_service;              /**< nC-array with (average) customer service times (in minutes; nonneg./pos. integers in [0,35]) */
   int**                 estim_t_service;        /**< (nC x 2)-array with estimated customer service time intervals (in minutes; pos. integers in [1,25]) */
   SCIP_Real*            driver_bond;            /**< nC-array with values specifying customers' sensitivity to change in driver (doubles in [0,1]) */
   int*                  t_walk;                 /**< nC-array with times driver needs to walk from truck to customers' house (in minutes; nonneg. integers in [0,5]) */
   int*                  grade;                  /**< nC-array with customer grade/value categories (1-8 representing strings A,B,C1,C2,C3,D1,D2,D3) */
   int*                  revenue;                /**< nC-array with customer revenue values (in EUR; nonneg./pos. integers in [0,5000]) */
   SCIP_Real*            riskfactor;             /**< nC-array with customer drop-out risk likelihood values (doubles in [0,1]) */
   int**                 distAM;                 /**< (nC x nC)-array with morning's distances between customers in data set (in meters; nonneg./pos. integers) */
   int**                 distNoon;               /**< (nC x nC)-array with noon's distances between customers in data set (in meters; nonneg./pos. integers) */
   int**                 distPM;                 /**< (nC x nC)-array with afternoon's distances between customers in data set (in meters; nonneg./pos. integers) */
   int**                 t_travelAM;             /**< (nC x nC)-array with morning's travel times between customers in data set (in seconds; nonneg./pos. integers) */
   int**                 t_travelNoon;           /**< (nC x nC)-array with noon's travel times between customers in data set (in seconds; nonneg./pos. integers) */
   int**                 t_travelPM;             /**< (nC x nC)-array with afternoon's travel times between customers in data set (in seconds; nonneg./pos. integers) */
   SCIP_Real*            depot_coordinates;      /**< 2-array with depot coordinates (doubles) */
   int*                  depot_out_distAM;       /**< nC-array with morning's distances from depot to customers in data set (in meters; nonneg./pos. integers) */
   int*                  depot_out_distNoon;     /**< nC-array with noon's distances from depot to customers in data set (in meters; nonneg./pos. integers) */
   int*                  depot_out_distPM;       /**< nC-array with afternoon's distances from depot to customers in data set (in meters; nonneg./pos. integers) */
   int*                  depot_out_t_travelAM;   /**< nC-array with morning's travel times from depot to customers in data set (in seconds; nonneg./pos. integers) */
   int*                  depot_out_t_travelNoon; /**< nC-array with noon's travel times from depot to customers in data set (in seconds; nonneg./pos. integers) */
   int*                  depot_out_t_travelPM;   /**< nC-array with afternoon's travel times from depot to customers in data set (in seconds; nonneg./pos. integers) */
   int*                  depot_in_distAM;        /**< nC-array with morning's distances from customers in data set to depot (in meters; nonneg./pos. integers) */
   int*                  depot_in_distNoon;      /**< nC-array with noon's distances from customers in data set to depot (in meters; nonneg./pos. integers) */
   int*                  depot_in_distPM;        /**< nC-array with afternoon's distances from customers in data set to depot (in meters; nonneg./pos. integers) */
   int*                  depot_in_t_travelAM;    /**< nC-array with morning's travel times from customers in data set to depot (in seconds; nonneg./pos. integers) */
   int*                  depot_in_t_travelNoon;  /**< nC-array with noon's travel times from customers in data set to depot (in seconds; nonneg./pos. integers) */
   int*                  depot_in_t_travelPM;    /**< nC-array with afternoon's travel times from customers in data set to depot (in seconds; nonneg./pos. integers) */
   SCIP_Bool             workOnSaturdays;        /**< indicator whether saturdays are to be counted as workdays or not (TRUE/FALSE (1/0 unsigned int.)) */
   date*                 startDate;              /**< struct containing the start date of the planning period (day, month, year, weekday) */
   date*                 endDate;                /**< struct containing the end date of the planning period (day, month, year, weekday) */
   date**                lastServiceDate;        /**< nC-array of structs containing the dates at which customers were last serviced (day, month, year, weekday) */
   dayProperty**         days;                   /**< struct containing the properties of each day */
   /* delivery dependencies (given customers must be serviced directly before or after current customer) -> ignored for now; precise meaning needs to be discussed further */
}instance_data;

/** model data (stores processed instance data, in the form that is used/required by the optimization model) */
typedef struct _modelData
{
   int                   nC;                     /**< number of customers + 1 (depot stored as last "customer") */
   int*                  customerIDs;            /**< nC-array with customer IDs (not needed for model, but kept for output/information purposes; depot gets ID -1) */
   modelWindow**         timeWindows;            /**< nC-array with linked list structures to store customers' delivery time windows */
   SCIP_Real*            obj;                    /**< nC-array with objective function coefficients for each customer */
   int                   shift_start;            /**< time when driver's daily work shift starts (no departure from depot before then) (in seconds) */
   int                   shift_end;              /**< time when driver's daily work shift ends (must return to depot by then) (in seconds) */
   int*                  t_service;              /**< nC-array with (average) customer service times (in seconds) */
   int*                  t_service_maxDev;       /**< nC-array with estimated maximal deviation of customer service times beyond average (in seconds) */
   int**                 t_travelAM;             /**< (nC x nC)-array with morning's travel times between customers in data set (in seconds) */
   int**                 t_travelNoon;           /**< (nC x nC)-array with noon's travel times between customers in data set (in seconds) */
   int**                 t_travelPM;             /**< (nC x nC)-array with afternoon's travel times between customers in data set (in seconds) */
   int**                 t_travel;               /**< (nC x nC)-array with average travel times between customers in data set (in seconds) */
   int**                 t_travel_maxDev;        /**< (nC x nC)-array with estimated maximal deviation of travel times between customers beyond average (in seconds) */
   SCIP_Bool             workOnSaturdays;        /**< indicator whether saturdays are to be counted as workdays or not (TRUE/FALSE (1/0 unsigned int.)); for information/output purposes */
   date*                 startDate;              /**< struct containing the start date of the planning period (day, month, year, weekday); for information/output purposes */
   date*                 endDate;                /**< struct containing the end date of the planning period (day, month, year, weekday); for information/output purposes */
   int                   maxDelayEvents;         /**< maximal number of delay events (deviations in service and travel times combined) to be "robustified" against */
   int                   nDays;                  /**< number of work days in planning period (from startDate to endDate, including both) */
   neighbor***           neighbors;              /**< (nC x nDays) array of linked lists with neighbor nodes of customer nodes in total network (can be NULL <-> complete digraph) */
   dayProperty**         days;                   /**< struct containing the properties of each day */
   int*                  day_sizes;              /**< nDays-array with the number of available customers for each day */
   /* delivery dependencies (given customers must be serviced directly before or after current customer) -> ignored for now; precise meaning needs to be discussed further */
}model_data;


/** help function for extracting values from json */
int get_long_from_json_obj(FILE** fp, char* scanbuf, const char* format, long* valuebuf, int to_read);

/** read solution json data file */
SCIP_RETCODE readSolutionData(SCIP* scip, const char* path, solution_data* solution);

/** free solution data structure */
SCIP_RETCODE deinitSolutionData(SCIP* scip, solution_data* solution);

SCIP_RETCODE writeSolutionJSON(SCIP* scip, solution_data* solutionData, char* outputFile);

/** read customer csv data file */
extern
SCIP_RETCODE readCustomerData(SCIP* scip, const char* customerFile, instance_data* instanceData);



/** read distance and travel-time data files */
extern
SCIP_RETCODE readTravelData(SCIP* scip, const char* matresponseAMFile, const char* matresponseNoonFile, const char* matresponsePMFile, instance_data* instanceData);



/** read meta data file */
extern
SCIP_RETCODE readMetaData(SCIP* scip, const char* metaFile, instance_data* instanceData);



/* set weekday for given date (by comparing with reference date macros) */
extern
SCIP_RETCODE setWeekday(int day,int month,int year,int* weekday);



/** initialize instance data struct */
extern
SCIP_RETCODE initInstanceData(SCIP* scip, instance_data* instanceData, int nC);



/** de-initialize instance data struct (clean up memory) */
extern
SCIP_RETCODE deinitInstanceData(SCIP* scip, instance_data* instanceData);



/** create (initialize and fill) model data struct */
extern
SCIP_RETCODE createModelData(SCIP* scip, model_data* modelData, instance_data* instanceData);



/** de-initialize model data struct (clean up memory) */
extern
SCIP_RETCODE deinitModelData(SCIP* scip, model_data* modelData);



/** derive and set model service times and assoc. max. deviations from instance data */
extern
SCIP_RETCODE setModelServiceTimes(int* t_service, int* t_service_maxDev, int twalk, int ts, int* ets);



/** derive and set model objective function coefficient for customer i from instance data */
extern
SCIP_RETCODE setModelObjectiveCoeff(SCIP_Real* obj, int i, instance_data* instanceData);



/** derive and set model delivery time windows for customer i from instance data */
extern
SCIP_RETCODE setModelTimeWindows(SCIP* scip, modelWindow** tw, int i, instance_data* instanceData);



/* increment the date in a date struct by a given number of days */
extern
SCIP_RETCODE incrementDate(date* d, int k);



/* decrement the date in a date struct by a given number of days */
extern
SCIP_RETCODE decrementDate(date* d, int k);



/* set model time window (work)day index */
extern
SCIP_RETCODE setModelDay(int* day, date* twDate, instance_data* instanceData);



/* get workday difference for given date pair */
extern
int getWorkdayDiff(date* date1, date* date2, instance_data* instanceData);



/* write model data to file */
extern
SCIP_RETCODE writeModelData(model_data* modelData, const char* outputFile);



/* read model data from file */
extern
SCIP_RETCODE readModelData(SCIP* scip, const char* inputFile, model_data* modelData);


#endif
