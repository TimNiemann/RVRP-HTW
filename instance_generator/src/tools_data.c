/**@file   tools_data.c
 * @brief  functions for handling of data (reading, writing, processing)
 * @author Andreas M. Tillmann, TU Braunschweig
 */

#include "tools_data.h"
#include <string.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <scip/scip.h>

/*#define SCIP_DEBUG*/ /* uncomment for some (debugging) output */

#define DISABLE_ALL_LASTSERVICEDATES FALSE /**< (TRUE/FALSE) toggle disabling last service date data (-> more flexible planning) */

const int SCAN_BUF_SIZE = 500;
const int VALUE_BUF_SIZE = 300;

/** utility function */
int get_long_from_json_obj(
        FILE** fp,
        char* scanbuf,
        const char* format,
        long* valuebuf,
        int to_read
        )
{
    char* iter;
    int count = 0, confirmed = 0;
    if(to_read == 0) return 0;
    while ((fgets(scanbuf, SCAN_BUF_SIZE * sizeof(char), *fp) != NULL)) {
        if (strstr(scanbuf, format) != NULL || confirmed == 1) {
            confirmed = 1;
            iter = scanbuf;
            while (*iter) {
                if (isdigit(*iter) || ((*iter == '-' || *iter == '+') && isdigit(*(iter + 1)))) {
                    valuebuf[count] = strtol(iter, &iter, 10);
                    count++;
                }
                else iter++;
            }
            if (count == to_read) return 0;
        }
    }
    printf("format error detected\n");
    return 1;
}

/** read solution json file */
SCIP_RETCODE readSolutionData(
        SCIP* scip,                     /** SCIP pointer */
        const char* path,               /** Input file name */
        solution_data* solution         /** pointer to solution data struct (alloc'd outside of this function!) */
        )
{
    long* valuebuf = malloc(VALUE_BUF_SIZE * sizeof(long));
    char* scanbuf = malloc(SCAN_BUF_SIZE * sizeof(char));
    FILE* fp = fopen(path, "r");
    if (fp == NULL) {
        printf("Unable to open %s\n", path);
        return 1;
    }
    if(get_long_from_json_obj(&fp ,scanbuf, "number of days", valuebuf, 1) != 0) return SCIP_READERROR;
    solution->nDays = (int)valuebuf[0];
    if(get_long_from_json_obj(&fp ,scanbuf, "number of customers", valuebuf, 1) != 0) return SCIP_READERROR;
    solution->nC = (int)valuebuf[0]+1;
    if(get_long_from_json_obj(&fp ,scanbuf, "robustness degree", valuebuf, 1) != 0) return SCIP_READERROR;
    solution->gamma = (int)valuebuf[0];
    if(get_long_from_json_obj(&fp ,scanbuf, "window length", valuebuf, 1) != 0) return SCIP_READERROR;
    solution->windowLength = (int)valuebuf[0];
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->alphas),3) );
    if(get_long_from_json_obj(&fp ,scanbuf, "alpha_one", valuebuf, 2) != 0) return SCIP_READERROR;
    solution->alphas[0] = valuebuf[0];
    if(valuebuf[1] != 0) solution->alphas[0] += (double)valuebuf[1] / pow(10, floor(log10((double)valuebuf[1]) + 1));
    if(get_long_from_json_obj(&fp ,scanbuf, "alpha_two", valuebuf, 2) != 0) return SCIP_READERROR;
    solution->alphas[1] = valuebuf[0];
    if(valuebuf[1] != 0) solution->alphas[1] += (double)valuebuf[1] / pow(10, floor(log10((double)valuebuf[1]) + 1));
    if(get_long_from_json_obj(&fp ,scanbuf, "alpha_three", valuebuf, 2) != 0) return SCIP_READERROR;
    solution->alphas[2] = valuebuf[0];
    if(valuebuf[1] != 0) solution->alphas[2] += (double)valuebuf[1] / pow(10, floor(log10((double)valuebuf[1]) + 1));
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->timeWindows),solution->nC) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->tourLength),solution->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->tourDuration),solution->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->customerIDs),solution->nC) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->obj),solution->nDays) );
    SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(solution->tour),solution->nDays) );
    for(int i=0; i<solution->nC - 1; i++) {
        SCIP_CALL( SCIPallocBlockMemory(scip, &(solution->timeWindows[i])) );
        SCIP_CALL( SCIPallocBlockMemory(scip, &(solution->timeWindows[i]->modelwindow)));
        if(get_long_from_json_obj(&fp ,scanbuf, "customerID", valuebuf, 9) != 0) return SCIP_READERROR;
        solution->customerIDs[i] = (int)valuebuf[0];
        solution->timeWindows[i]->day = (int)valuebuf[1];
        solution->timeWindows[i]->start_t = (int)valuebuf[2];
        solution->timeWindows[i]->end_t = (int)valuebuf[3];
        solution->timeWindows[i]->latest_arrival = (int)valuebuf[4];
        solution->timeWindows[i]->modelwindow->next = NULL;
        solution->timeWindows[i]->modelwindow->day = (int)valuebuf[1];
        solution->timeWindows[i]->modelwindow->start_t = (int)valuebuf[5];
        solution->timeWindows[i]->modelwindow->end_t = (int)valuebuf[6];
        solution->timeWindows[i]->modelwindow->weigth = (double)valuebuf[7];
        if(valuebuf[8] != 0) solution->timeWindows[i]->modelwindow->weigth += (double)valuebuf[8] / pow(10, floor(log10((double)valuebuf[8]) + 1));
    }
    solution->customerIDs[solution->nC - 1] = -1;
    for(int i=0; i<solution->nDays; i++){
        if(get_long_from_json_obj(&fp ,scanbuf, "day", valuebuf, 5) != 0) return SCIP_READERROR;
        solution->tourLength[i] = (int)valuebuf[1];
        solution->tourDuration[i] = (int)valuebuf[2];
        solution->obj[i] = (double)valuebuf[3];
        if(valuebuf[4] != 0) solution->obj[i] += (double)valuebuf[4] / pow(10, floor(log10((double)valuebuf[4]) + 1));
        if(get_long_from_json_obj(&fp ,scanbuf, "visited customers", valuebuf, solution->tourLength[i]) != 0) return SCIP_READERROR;
        SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(solution->tour[i]), solution->tourLength[i]) );
        for(int j=0; j<solution->tourLength[i]; j++)
            solution->tour[i][j] = (int)valuebuf[j];
    }
    fclose(fp);
    free(scanbuf);
    free(valuebuf);
    return SCIP_OKAY;
}

SCIP_RETCODE deinitSolutionData(
        SCIP* scip,
        solution_data* solution
        )
{
    for (int i=(solution->nDays)-1; i >= 0; i--){
        SCIPfreeBlockMemoryArray(scip, &(solution->tour[i]), solution->tourLength[i]);
    }
    for (int i=(solution->nC)-2; i >= 0; i--){
        SCIPfreeBlockMemory(scip, &(solution->timeWindows[i]->modelwindow));
        SCIPfreeBlockMemory(scip, &(solution->timeWindows[i]));
    }
    SCIPfreeBlockMemoryArray(scip, &(solution->alphas), 3);
    SCIPfreeBlockMemoryArray(scip, &(solution->timeWindows), solution->nC);
    SCIPfreeBlockMemoryArray(scip, &(solution->tour), solution->nDays);
    SCIPfreeBlockMemoryArray(scip, &(solution->tourLength), solution->nDays);
    SCIPfreeBlockMemoryArray(scip, &(solution->tourDuration), solution->nDays);
    SCIPfreeBlockMemoryArray(scip, &(solution->customerIDs), solution->nC);
    SCIPfreeBlockMemoryArray(scip, &(solution->obj), solution->nDays);
    return SCIP_OKAY;
}


SCIP_RETCODE writeSolutionJSON(
        SCIP*               scip,
        solution_data*      solutionData,
        char*               outputFile
)
{
    FILE *fp;
    int i;

    assert(solutionData != NULL);
    assert(outputFile != NULL);

    fp = fopen(outputFile, "w");
    if (fp == NULL)
    {
        SCIPwarningMessage(scip, "Can't open outputfile: %s\n", outputFile);
        return SCIP_WRITEERROR;
    }

    fprintf(fp, "{\n");
    fprintf(fp, "\t\"number of days\": %d,\n", solutionData->nDays);
    fprintf(fp, "\t\"number of customers\": %d,\n", solutionData->nC - 1);
    fprintf(fp, "\t\"robustness degree\": %d,\n", solutionData->gamma);
    fprintf(fp, "\t\"window length\": %d,\n", solutionData->windowLength);
    fprintf(fp, "\t\"alpha_one\": %.2f,\n", solutionData->alphas[0]);
    fprintf(fp, "\t\"alpha_two\": %.2f,\n", solutionData->alphas[1]);
    fprintf(fp, "\t\"alpha_three\": %.2f,\n", solutionData->alphas[2]);
    fprintf(fp, "\t\"customers\": [\n");
    for (i = 0; i < solutionData->nC - 1; i++)
    {
        fprintf(fp, "\t\t{\n");
        fprintf(fp, "\t\t\t\"customerID\": %d,\n", solutionData->customerIDs[i]);
        fprintf(fp, "\t\t\t\"delivery day\": %d,\n", solutionData->timeWindows[i]->day);
        fprintf(fp, "\t\t\t\"delivery window start\": %d,\n", solutionData->timeWindows[i]->start_t);
        fprintf(fp, "\t\t\t\"delivery window end\": %d,\n", solutionData->timeWindows[i]->end_t);
        fprintf(fp, "\t\t\t\"delivery latest arrival\": %d,\n", solutionData->timeWindows[i]->latest_arrival);
        fprintf(fp, "\t\t\t\"hard window start\": %d,\n", solutionData->timeWindows[i]->modelwindow->start_t);
        fprintf(fp, "\t\t\t\"hard window end\": %d,\n", solutionData->timeWindows[i]->modelwindow->end_t);
        fprintf(fp, "\t\t\t\"encounter probability\": %.2f\n", solutionData->timeWindows[i]->modelwindow->weigth);
        fprintf(fp, "\t\t}");
        if (i != solutionData->nC - 2)
        {
            fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\t],\n");

    fprintf(fp, "\t\"tours\": [\n");
    for (i = 0; i < solutionData->nDays; i++)
    {
        int j;
        fprintf(fp, "\t\t{\n");
        fprintf(fp, "\t\t\t\"day\": %d,\n", i);
        fprintf(fp, "\t\t\t\"length\": %d,\n", solutionData->tourLength[i]);
        fprintf(fp, "\t\t\t\"duration\": %d,\n", solutionData->tourDuration[i]);
        fprintf(fp, "\t\t\t\"objective\": %f,\n", solutionData->obj[i]);
        fprintf(fp, "\t\t\t\"visited customers\": [\n");
        for (j = 0; j < solutionData->tourLength[i]; j++)
        {
            fprintf(fp, "\t\t\t\t%d", solutionData->tour[i][j]);
            if (j != solutionData->tourLength[i] - 1)
            {
                fprintf(fp, ",");
            }
            fprintf(fp, "\n");
        }
        fprintf(fp, "\t\t\t]\n");
        fprintf(fp, "\t\t}");
        if (i != solutionData->nDays - 1)
        {
            fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fprintf(fp, "\t]\n");
    fprintf(fp, "}\n");

    fclose(fp);
    printf("\nOutputfile created: %s\n", outputFile);
    return SCIP_OKAY;
}


/* read customers.csv file */
SCIP_RETCODE readCustomerData(
   SCIP* scip,                            /**< SCIP pointer */
   const char* customerFile,              /**< customer data file name (full path) */
   instance_data* instanceData            /**< struct in which to store the data read from file */
   )
{
   FILE* inFILE = NULL;
   int i,j,k,l;
   char strbuffer[SCIP_MAXSTRLEN];
   char strbuffer2[SCIP_MAXSTRLEN];
   int success;
   int itmp;
   timeWindow* tmp_tw = NULL;
   timeWindow* twHead = NULL;

   assert( scip != NULL );
   assert( customerFile != NULL );
   assert( instanceData != NULL );

   inFILE = fopen(customerFile,"r");
   assert( inFILE != NULL );

   /* get data set size (number of lines in csv file) */
   k = -1; /* init. with -1 so that header line is not counted */
   while( fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      k++;
   }
   assert( k > 0 );
   rewind(inFILE);

   /* now that k, the data set size, is known, memory can be allocated accordingly: */
   SCIP_CALL( initInstanceData(scip, instanceData, k) );

   /* now, scan the file again and extract and store the desired data: */
   /* NB: customer data csv files are assumed to contain the following 28 columns (in that order!); columns with names in
      square brackets will be skipped when reading the file as they are (currently) not needed for the optimization */
   /*
     Kundennummer; [Straße]; [Hausnummer]; [PLZ]; [Ort]; Verortung; [Lieferwunschstartzeit]; Zeitfenster,
     [Lieferwunschende]; Belieferungsrhythmus, Erreichbarkeit, Erreichbarkeitsvermutung, Terminsensibilität,
     Termintoleranz, [verschobener Termin]; [Nachlieferung]; Bearbeitungszeit, Verweildauervermutung,
     Verkaufsfahrerbindung, Wegzeit Fzg, Kundenwertigkeit, Kundenumsatzwert, [Warenpräferenz]; Risikofaktor,
     [Kampagnenkunde]; [Neukunde]; Lieferabhängigkeit; Letzte Lieferung
   */
   /* Lieferabhängigkeit is currently still being skipped -> we can change that in the future, if so desired */
   /* first line (headers) -> ignore/skip */
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading first line.\n");
      return SCIP_READERROR;
   }

   /* read data from remaining lines in file: */
   l = 0; /* to count number of data lines read */
   while( fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      k = 0; /* to keep track of current position in string containing whole line */
      j = 0; /* to count number of processed fields (columns) -> check that this is 27 in each processed line */
      success = 0; /* to detect empty fields */
      while( k < (int)strlen(strbuffer) )
      {
         success = sscanf(&strbuffer[k], "%[^;]", strbuffer2); /* extract field content up to next semicolon */
         j++;
         if( success ) /* field was not empty -> process contents */
         {
            k += (int)strlen(strbuffer2) +1;
            switch(j)
            {
            case 1: /* Kundennummer; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->customerIDs[l] = atoi(strbuffer2);
               break;
            case 6: /* Verortung; assumed to be given as a string of the form POINT(a b), with a and b floating point numbers (with 6 decimals) */
               itmp = sscanf(&strbuffer2[6], "%lf %lf", &instanceData->coordinates[l][0], &instanceData->coordinates[l][1]);
               assert( itmp == 2 );
               break;
            case 8: /* Zeitfenster; assumed to be provided in one of the following two formats:
                          ['Mo 15:00-19:00', 'Mi 15:00-18:00', 'Fr 16:00-17:00']     (for example), or
                          06.10.2020 - 07.10.2020 14:00-15:00
                       The first case is interpreted as time windows for each specified weekday over the whole
                       planning period (start date to end date), subject to delivery frequency w.r.t. last
                       service date. (German weekday abbreviations 'Mo'-'Sa' <-> 'Mon'-'Sat' <-> 0-6.)
                       The second case is interpreted as time windows for each day over the specified subperiod
                       of the whole planning period. NB: In the 2nd case, delivery frequency is ignored (or
                       equivalently, assumed to be 0).
                       NB: Non-compliance with the above formatting restrictions cannot be handled -> error.
                     */
               i = 0;
               twHead = instanceData->timeWindows[l]; /* (init.) pointer to current list element */
               tmp_tw = NULL;
               assert( twHead == NULL ); /* first list element was initialized with NULL, to allow for entirely empty lists */
               while( (strncmp(&strbuffer2[i],"[\'",2) == 0) || (strncmp(&strbuffer2[i],"\',",2) == 0) ) /* first case: given list of weekday-bound time windows */
               {
                  i += 2;
                  if( strncmp(&strbuffer2[i]," \'",2) == 0 ) /* space between list entries -> skip */
                  {
                     i += 2;
                  }
                  /* allocate memory for new list element, then parse weekday and time window data point: */
                  SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
                  assert( tmp_tw != NULL );
                  tmp_tw->recurring = TRUE;
                  SCIP_CALL( SCIPallocBlockMemory(scip, &(tmp_tw->twDate)) );
                  assert( tmp_tw->twDate != NULL );
                  (tmp_tw->twDate)->day   = -1; /* set d/m/y to -1, since only weekday matters */
                  (tmp_tw->twDate)->month = -1;
                  (tmp_tw->twDate)->year  = -1;
                  if( strncmp(&strbuffer2[i],"Mo ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 0;
                  }
                  else if( strncmp(&strbuffer2[i],"Di ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 1;
                  }
                  else if( strncmp(&strbuffer2[i],"Mi ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 2;
                  }
                  else if( strncmp(&strbuffer2[i],"Do ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 3;
                  }
                  else if( strncmp(&strbuffer2[i],"Fr ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 4;
                  }
                  else if( strncmp(&strbuffer2[i],"Sa ",3) == 0 )
                  {
                     (tmp_tw->twDate)->weekday = 5;
                  }
                  else
                  {
                     itmp = sscanf(&strbuffer2[i], "%d.%d.%d",&(tmp_tw->twDate->day),&(tmp_tw->twDate->month),&(tmp_tw->twDate->year));
                     if(itmp != 3)
                     {
                         fprintf(stderr,"could not read/interpret Zeitfenster in line %d\n",l+1);
                         return SCIP_READERROR;
                     }
                     SCIP_CALL( setWeekday(tmp_tw->twDate->day,tmp_tw->twDate->month,tmp_tw->twDate->year,&(tmp_tw->twDate->weekday)) );
                     tmp_tw->recurring = FALSE;
                     i += 8;
                  }
                  i += 3;
                  itmp = sscanf(&strbuffer2[i], "%d:%d-%d:%d",&(tmp_tw->start_hour),&(tmp_tw->start_min),&(tmp_tw->end_hour),&(tmp_tw->end_min));
                  assert( itmp == 4 );
                  i += 9; /* (h)h:mm-(h)h:mm has 9,10, or 11 characters -> move position accordingly */
                  while( strncmp(&strbuffer2[i], "\'", 1) != 0 && strncmp(&strbuffer2[i], "|", 1) != 0 )
                  {
                     i++;
                  }
                  if(strncmp(&strbuffer2[i], "|", 1) == 0) /* save weight if given*/
                  {
                     i++;
                     itmp = sscanf(&strbuffer2[i], "%lf", &(tmp_tw->weigth));
                     assert(itmp == 1);
                     while( strncmp(&strbuffer2[i], "\'", 1) != 0)
                     {
                        i++;
                     }
                  }else /* default weight = 1.0 */
                  {
                     tmp_tw->weigth = 1.0;
                  }
                  tmp_tw->next = NULL;
                  if( twHead == NULL ) /* list element just constructed is the first one */
                  {
                     twHead = tmp_tw;
                     instanceData->timeWindows[l] = tmp_tw;
                     tmp_tw = tmp_tw->next;
                     assert( tmp_tw == NULL );
                     assert( instanceData->timeWindows[l] != NULL );
                  }
                  else /* list element just constructed gets appended to end of the current list */
                  {
                     twHead->next = tmp_tw;
                     twHead = twHead->next;
                     tmp_tw = tmp_tw->next;
                     assert( tmp_tw == NULL );
                     assert( twHead != NULL );
                  }
               }
               if( twHead != NULL )
               {
                  assert( i > 0 );
                  assert( strncmp(&strbuffer2[i],"\']",2) == 0 ); /* must have reached end of given list */
               }
               else /* second case: given concrete date range with time windows applicable to each day thereof */
               {
                  /* in this case, we know the exact format of the string to process, and know that we will need just
                     two list elements, so we allocate and fill them directly: */
                  SCIP_CALL( SCIPallocBlockMemory(scip, &(instanceData->timeWindows[l])) );
                  assert( instanceData->timeWindows[l] != NULL );
                  SCIP_CALL( SCIPallocBlockMemory(scip, &((instanceData->timeWindows[l])->next)) );
                  assert( (instanceData->timeWindows[l])->next != NULL );
                  (instanceData->timeWindows[l])->recurring = FALSE;
                  ((instanceData->timeWindows[l])->next)->recurring = FALSE;
                   (instanceData->timeWindows[l])->weigth = 1.0;
                   (instanceData->timeWindows[l])->next->weigth = 1.0;
                  SCIP_CALL( SCIPallocBlockMemory(scip, &((instanceData->timeWindows[l])->twDate)) );
                  assert( (instanceData->timeWindows[l])->twDate != NULL );
                  SCIP_CALL( SCIPallocBlockMemory(scip, &(((instanceData->timeWindows[l])->next)->twDate)) );
                  assert( ((instanceData->timeWindows[l])->next)->twDate != NULL );
                  ((instanceData->timeWindows[l])->next)->next = NULL;
                  /* format is like this: "06.10.2020 - 07.10.2020 14:00-15:00" (without quotation marks) */
                  itmp = sscanf(strbuffer2, "%d.%d.%d - %d.%d.%d %d:%d-%d:%d", &(((instanceData->timeWindows[l])->twDate)->day), &(((instanceData->timeWindows[l])->twDate)->month), &(((instanceData->timeWindows[l])->twDate)->year), &((((instanceData->timeWindows[l])->next)->twDate)->day), &((((instanceData->timeWindows[l])->next)->twDate)->month), &((((instanceData->timeWindows[l])->next)->twDate)->year), &((instanceData->timeWindows[l])->start_hour), &((instanceData->timeWindows[l])->start_min), &((instanceData->timeWindows[l])->end_hour),  &((instanceData->timeWindows[l])->end_min));
                  assert( itmp == 10 );
                  SCIP_CALL( setWeekday(((instanceData->timeWindows[l])->twDate)->day,((instanceData->timeWindows[l])->twDate)->month,((instanceData->timeWindows[l])->twDate)->year,&(((instanceData->timeWindows[l])->twDate)->weekday)) );
                  SCIP_CALL( setWeekday((((instanceData->timeWindows[l])->next)->twDate)->day,(((instanceData->timeWindows[l])->next)->twDate)->month,(((instanceData->timeWindows[l])->next)->twDate)->year,&((((instanceData->timeWindows[l])->next)->twDate)->weekday)) );
                  /* time window start and end times are the same in second element as in the first one */
                  ((instanceData->timeWindows[l])->next)->start_hour = (instanceData->timeWindows[l])->start_hour;
                  ((instanceData->timeWindows[l])->next)->start_min  = (instanceData->timeWindows[l])->start_min;
                  ((instanceData->timeWindows[l])->next)->end_hour   = (instanceData->timeWindows[l])->end_hour;
                  ((instanceData->timeWindows[l])->next)->end_min    = (instanceData->timeWindows[l])->start_min;
               }
               break;
            case 10: /* Belieferungsrhythmus; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->deliveryFreq[l] = atoi(strbuffer2);
               break;
            case 11: /* Erreichbarkeit; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->availability[l] = atoi(strbuffer2);
               break;
            case 12: /* Erreichbarkeitsvermutung; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->estim_availability[l] = (SCIP_Real)atof(strbuffer2);
               break;
            case 13: /* Terminsensibilität; assumed to be given directly as a (single) character, without quotation marks or other characters */
               if( strcmp(strbuffer2,"A")==0 )
               {
                  instanceData->appointment_sens[l] = 1;
               }
               else if( strcmp(strbuffer2,"B")==0 )
               {
                  instanceData->appointment_sens[l] = 2;
               }
               else if( strcmp(strbuffer2,"C")==0 )
               {
                  instanceData->appointment_sens[l] = 3;
               }
               else
               {
                  fprintf(stderr,"could not read/interpret Terminsensibilitaet in line %d\n",l+1);
                  return SCIP_READERROR;
               }
               break;
            case 14: /* Termintoleranz; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->appointment_tol[l] = atoi(strbuffer2);
               break;
            case 17: /* Bearbeitungszeit; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->t_service[l] = atoi(strbuffer2);
               break;
            case 18: /* Verweildauervermutung; assumed to be given as a string of the form [a, b], with a and b integer numbers */
               itmp = sscanf(&strbuffer2[1], "%d, %d", &instanceData->estim_t_service[l][0], &instanceData->estim_t_service[l][1]);
               assert( itmp == 2 );
               break;
            case 19: /* Verkaufsfahrerbindung; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->driver_bond[l] = (SCIP_Real)atof(strbuffer2);
               break;
            case 20: /* Wegzeit Fzg; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->t_walk[l] = atoi(strbuffer2);
               break;
            case 21: /* Kundenwertigkeit; assumed to be given as a string consisting of 1 or 2 characters, without quotation marks or other characters */
               if( strcmp(strbuffer2,"A")==0 )
               {
                  instanceData->grade[l] = 1;
               }
               else if( strcmp(strbuffer2,"B")==0 )
               {
                  instanceData->grade[l] = 2;
               }
               else if( strncmp(strbuffer2,"C",1)==0 )
               {
                  instanceData->grade[l] = 2+atoi(&strbuffer2[1]); /* C1->3, C2->4, C3->5 */
               }
               else if (strncmp(strbuffer2,"D",1)==0 )
               {
                  instanceData->grade[l] = 5+atoi(&strbuffer2[1]); /* D1->6, D2->7, D3->8 */
               }
               else
               {
                  fprintf(stderr,"could not read/interpret Kundenwertigkeit in line %d\n",l+1);
                  return SCIP_READERROR;
               }
               break;
            case 22: /* Kundenumsatzwert; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->revenue[l] = atoi(strbuffer2);
               break;
            case 24: /* Risikofaktor; assumed to be given directly as a number without quotation marks or other characters */
               instanceData->riskfactor[l] = (SCIP_Real)atof(strbuffer2);
               break;
            /* case 27: Lieferabhängigkeit; ... currently skipped */
            case 28: /* Letzte Lieferung; date of last service at customer */
               /* check if empty (or whether use of last service dates has been disabled): */
               if( (strncmp(strbuffer2, "[]", 2) == 0) || (DISABLE_ALL_LASTSERVICEDATES) )
                  {
                  (instanceData->lastServiceDate[l])->day     = -1;
                  (instanceData->lastServiceDate[l])->month   = -1;
                  (instanceData->lastServiceDate[l])->year    = -1;
                  (instanceData->lastServiceDate[l])->weekday = -1;
               }
               else
               {
                  itmp = sscanf(strbuffer2, "%d.%d.%d", &((instanceData->lastServiceDate[l])->day),&((instanceData->lastServiceDate[l])->month),&((instanceData->lastServiceDate[l])->year));
                  assert( itmp == 3 );
                  SCIP_CALL( setWeekday((instanceData->lastServiceDate[l])->day,(instanceData->lastServiceDate[l])->month,(instanceData->lastServiceDate[l])->year,&((instanceData->lastServiceDate[l])->weekday)) );
               }
               break;
            default: /* fields 2,3,4,5,7,(8 for now),9,15,16,23,25,26,(27 for now) are skipped */
               break;
            }
         }
         else /* field was empty -> skip */
         {
            k++;
            /* possible TO DO: if column needs not contain values for every line, use defaults if empty field is encountered */
         }
      }
      l++;
      assert( j == 28 );
   }
   assert( l == instanceData->nC ); /* all data lines should have been processed */

   fclose(inFILE);

   return SCIP_OKAY;
}



/* read matrix response data files (matrix-response-{morning,noon,afternoon}.json) */
SCIP_RETCODE readTravelData(
   SCIP* scip,                            /**< SCIP pointer */
   const char* matresponseAMFile,         /**< morning distance & travel time (customer matrix) data file (full path) */
   const char* matresponseNoonFile,       /**< noon distance & travel time (customer matrix) data file (full path) */
   const char* matresponsePMFile,         /**< afternoon distance & travel time (customer matrix) data file (full path) */
   instance_data* instanceData             /**< struct in which to store the data read from file */
   )
{
   FILE* inFILE = NULL;
   int i,j,k,l;
   char strbuffer[SCIP_MAXSTRLEN];
   int success = 0;

   assert( scip != NULL );
   assert( matresponseAMFile != NULL );
   assert( matresponseNoonFile != NULL );
   assert( matresponsePMFile != NULL );
   assert( instanceData != NULL );

   /* NB: the file format is assumed to be as follows:
      the first 10 lines are syntax and meta-data and can hence be skipped,
      then, there will be nC (number of customers) many blocks consisting of
         line " "distance": [" -> skip
         nC distance values, one per line
         lines " ], " and " "duration: [" -> skip
         nC travel time values, one per line
         three lines of syntax -> skip
      the final line is again just syntax and can be skipped */
   /* THIS SYNTAX IS  N O T  CHECKED, but directly exploited for parsing */

   /* read morning distance & travel time data */
   inFILE = fopen(matresponseAMFile,"r");
   assert( inFILE != NULL );
   /* ignore/skip first 11 lines: */
   l = 0;
   while( l < 11 )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* traverse file and read distance and travel time data, skipping syntax-only lines */
   for (i = 0; i < instanceData->nC; i++) /* matrix rows */
   {
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: distance from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->distAM[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 2; k++) /* skip two syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: travel times from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->t_travelAM[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 4; k++) /* skip four syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   assert( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) ); /* EOF */
   fclose(inFILE);

   /* read noon distance & travel time data */
   inFILE = fopen(matresponseNoonFile,"r");
   assert( inFILE != NULL );
   /* ignore/skip first 11 lines: */
   l = 0;
   while( l < 11 )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* traverse file and read distance and travel time data, skipping syntax-only lines */
   for (i = 0; i < instanceData->nC; i++) /* matrix rows */
   {
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: distance from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->distNoon[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 2; k++) /* skip two syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: travel times from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->t_travelNoon[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 4; k++) /* skip four syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   assert( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) ); /* EOF */
   fclose(inFILE);

   /* read afternoon distance & travel time data */
   inFILE = fopen(matresponsePMFile,"r");
   assert( inFILE != NULL );
   /* ignore/skip first 11 lines: */
   l = 0;
   while( l < 11 )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* traverse file and read distance and travel time data, skipping syntax-only lines */
   for (i = 0; i < instanceData->nC; i++) /* matrix rows */
   {
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: distance from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->distPM[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 2; k++) /* skip two syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      for (j = 0; j < instanceData->nC; j++) /* columns */
      {
         l++;
         /* read nC values into matrix: travel times from customer i (row) to customer j (column) */
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error reading line %d.\n",l);
            return SCIP_READERROR;
         }
         success = sscanf(strbuffer, "%d", &(instanceData->t_travelPM[i][j]));
         assert( success == 1 );
      }
      for (k = 0; k < 4; k++) /* skip four syntax lines: */
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   assert( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) ); /* EOF */
   fclose(inFILE);

   return SCIP_OKAY;
}



/* read meta.json file */
SCIP_RETCODE readMetaData(
   SCIP* scip,                            /**< SCIP pointer */
   const char* metaFile,                  /**< meta data file (full path) */
   instance_data* instanceData             /**< struct in which to store the data read from file */
   )
{
   FILE* inFILE = NULL;
   int i,k,l;
   char strbuffer[SCIP_MAXSTRLEN];
   int success = 0;

   assert( scip != NULL );
   assert( metaFile != NULL );
   assert( instanceData != NULL );

   /* NB: the file format is assumed to be as follows:

      1 line syntax -> skip
          "startDate": "< DD.MM.YYYY >",
          "endDate": "< DD.MM.YYYY >",
      1 line syntax -> skip
              < double value: depot coordinate 1 >,
              < double value: depot coordinate 2 >
      1 line syntax -> skip
          "workOnSaturday": < true/false value (string) >,
      1 line syntax -> skip
      then, four blocks containing default/morning/noon/afternoon (in that order)
          distance and travel time values from depot to customers in the format:
            4 lines syntax -> skip
            nC distance values (from depot to customers)
            2 lines syntax -> skip
            nC travel time values (from depot to customers)
            4 lines syntax -> skip
      2 lines syntax -> skip
      then, four blocks containing default/morning/noon/afternoon (in that order)
          distance and travel time values from customers to depot in the format:
            2 lines syntax -> skip
            nC times: 2 lines syntax (->skip), distance value, 2 lines syntax (->skip), travel time value, 2 lines syntax (->skip)
            2 lines syntax -> skip
      the final 2 lines area again syntax and can be skipped */

   /* THIS SYNTAX IS  N O T  CHECKED, but directly exploited for parsing */

   /* read meta data */
   inFILE = fopen(metaFile,"r");
   assert( inFILE != NULL );
   /* ignore/skip first line: */
   l = 0;
   while( l < 1 )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* read startDate: */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading line %d.\n",l);
      return SCIP_READERROR;
   }
   /* find position of colon (from "startDate":) in strbuffer: */
   k = 0;
   while( strbuffer[k] != ':' )
   {
      k++;
   }
   assert( strbuffer[k] == ':' );
   success = sscanf(&strbuffer[k], ": \"%d.%d.%d", &((instanceData->startDate)->day), &((instanceData->startDate)->month), &((instanceData->startDate)->year));
   assert( success == 3 );
   SCIP_CALL( setWeekday((instanceData->startDate)->day,(instanceData->startDate)->month,(instanceData->startDate)->year,&((instanceData->startDate)->weekday)) );
   /* read endDate: */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading line %d.\n",l);
      return SCIP_READERROR;
   }
   /* find position of colon (from "endDate":) in strbuffer: */
   k = 0;
   while( strbuffer[k] != ':' )
   {
      k++;
   }
   assert( strbuffer[k] == ':' );
   success = sscanf(&strbuffer[k], ": \"%d.%d.%d", &((instanceData->endDate)->day), &((instanceData->endDate)->month), &((instanceData->endDate)->year));
   assert( success == 3 );
   SCIP_CALL( setWeekday((instanceData->endDate)->day,(instanceData->endDate)->month,(instanceData->endDate)->year,&((instanceData->endDate)->weekday)) );

   /* skip one syntax line */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error skipping line %d.\n",l);
      return SCIP_READERROR;
   }
   /* read first depot coordinate */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading line %d.\n",l);
      return SCIP_READERROR;
   }
   success = sscanf(strbuffer, "%lf", &(instanceData->depot_coordinates[0]));
   assert( success == 1 );
   /* read second depot coordinate */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading line %d.\n",l);
      return SCIP_READERROR;
   }
   success = sscanf(strbuffer, "%lf", &(instanceData->depot_coordinates[1]));
   assert( success == 1 );
   /* skip one syntax line */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error skipping line %d.\n",l);
      return SCIP_READERROR;
   }
   /* read workOnSaturdays value: */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error reading line %d.\n",l);
      return SCIP_READERROR;
   }
   /* find position of colon (from "workOnSaturday":) in strbuffer: */
   k = 0;
   while( strbuffer[k] != ':' )
   {
      k++;
   }
   assert( strbuffer[k] == ':' );
   if( strncmp(&strbuffer[k+2],"false",5) == 0 )
   {
      instanceData->workOnSaturdays = FALSE;
   }
   else
   {
      assert( strncmp(&strbuffer[k+2],"true",4) == 0 );
      instanceData->workOnSaturdays = TRUE;
   }
   /* skip one syntax line ("fromDepot") */
   l++;
   if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
   {
      fprintf(stderr,"Error skipping line %d.\n",l);
      return SCIP_READERROR;
   }
   /* skip 10+2*nC lines (syntax and "default" distances/travel times from depot to customers */
   for (k = 0; k < 10+2*(instanceData->nC); k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with morning distances and travel times from depot to customers */
   /* first, skip 4 syntax lines */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC distance values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_distAM[i]));
      assert( success == 1 );
   }
   /* skip two more syntax lines: */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC travel time values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_t_travelAM[i]));
      assert( success == 1 );
   }
   /* skip four more syntax lines, ending the block: */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with noon distances and travel times from depot to customers */
   /* first, skip 4 syntax lines */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC distance values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_distNoon[i]));
      assert( success == 1 );
   }
   /* skip two more syntax lines: */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC travel time values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_t_travelNoon[i]));
      assert( success == 1 );
   }
   /* skip four more syntax lines, ending the block: */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with afternoon distances and travel times from depot to customers */
   /* first, skip 4 syntax lines */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC distance values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_distPM[i]));
      assert( success == 1 );
   }
   /* skip two more syntax lines: */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now read nC travel time values */
   for (i = 0; i < instanceData->nC; i++ )
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_out_t_travelPM[i]));
      assert( success == 1 );
   }
   /* skip four more syntax lines, ending the block: */
   for (k = 0; k < 4; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* skip two syntax lines, ending the "fromDepot" and beginning the "toDepot" block: */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* skip 4+8*nC lines (syntax and "default" distances/travel times to depot from customers */
   for (k = 0; k < 4+8*(instanceData->nC); k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with morning distances and travel times to depot from customers */
   /* first, skip 2 syntax lines */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now, nC times, skip 2 lines, read distance value, skip 2 lines, read travel time value, skip 2 lines: */
   for (i = 0; i < instanceData->nC; i++)
   {
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_distAM[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_t_travelAM[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   /* finally, skip two syntax lines, ending the block */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with noon distances and travel times to depot from customers */
   /* first, skip 2 syntax lines */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now, nC times, skip 2 lines, read distance value, skip 2 lines, read travel time value, skip 2 lines: */
   for (i = 0; i < instanceData->nC; i++)
   {
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_distNoon[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_t_travelNoon[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   /* finally, skip two syntax lines, ending the block */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* parse block with afternoon distances and travel times to depot from customers */
   /* first, skip 2 syntax lines */
   for (k = 0; k < 2; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
   /* now, nC times, skip 2 lines, read distance value, skip 2 lines, read travel time value, skip 2 lines: */
   for (i = 0; i < instanceData->nC; i++)
   {
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_distPM[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error reading line %d.\n",l);
         return SCIP_READERROR;
      }
      success = sscanf(strbuffer, "%d", &(instanceData->depot_in_t_travelPM[i]));
      assert( success == 1 );
      for (k = 0; k < 2; k++)
      {
         l++;
         if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
         {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
         }
      }
   }
   /* finally, skip two syntax lines, ending the block */
   for (k = 0; k < 3; k++)
   {
      l++;
      if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
      {
         fprintf(stderr,"Error skipping line %d.\n",l);
         return SCIP_READERROR;
      }
   }
    l++;
    if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
    {
        fprintf(stderr,"Error skipping line %d.\n",l);
        return SCIP_READERROR;
    }
    l++;
    /* set nDays */
    SCIP_CALL( setModelDay(&(instanceData->nDays),instanceData->endDate,instanceData) ); /* model day-index of end date of planning period */
    instanceData->nDays += 1; /* now equals number of work days from start to end of planning period, including start and end date */
    if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) ) // no pause properties given
    {
        fclose(inFILE);
        instanceData->days = NULL;
    }
    else
    {
        int hour, minute;
        /* allocate space for day properties */
        SCIP_CALL( SCIPallocBlockMemoryArray( scip, &(instanceData->days), instanceData->nDays) );
        for(i = 0; i < instanceData->nDays; i++)
        {
            SCIP_CALL( SCIPallocBlockMemory( scip, &(instanceData->days[i])) );
            /* read index */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            instanceData->days[i]->index = i;
            /* read shift start */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            k = 0;
            while( strbuffer[k] != ':')
            {
                k++;
            }
            assert(strbuffer[k] == ':');
            success = sscanf(&strbuffer[k], ": \"%d:%d", &hour, &minute);
            assert(success == 2);
            instanceData->days[i]->shift_start = hour * 3600 + minute * 60;
            /* read shift end */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            k = 0;
            while( strbuffer[k] != ':')
            {
                k++;
            }
            assert(strbuffer[k] == ':');
            success = sscanf(&strbuffer[k], ": \"%d:%d", &hour, &minute);
            assert(success == 2);
            instanceData->days[i]->shift_end = hour * 3600 + minute * 60;
            /* read break start */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            k = 0;
            while( strbuffer[k] != ':')
            {
                k++;
            }
            assert(strbuffer[k] == ':');
            success = sscanf(&strbuffer[k], ": \"%d:%d", &hour, &minute);
            assert(success == 2);
            instanceData->days[i]->break_start = hour * 3600 + minute * 60;
            /* read break end */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            k = 0;
            while( strbuffer[k] != ':')
            {
                k++;
            }
            assert(strbuffer[k] == ':');
            success = sscanf(&strbuffer[k], ": \"%d:%d", &hour, &minute);
            assert(success == 2);
            instanceData->days[i]->break_end = hour * 3600 + minute * 60;
            /* read break length */
            l++;
            if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
            {
                fprintf(stderr,"Error skipping line %d.\n",l);
                return SCIP_READERROR;
            }
            k = 0;
            while( strbuffer[k] != ':')
            {
                k++;
            }
            assert(strbuffer[k] == ':');
            success = sscanf(&strbuffer[k], ": %d",  &minute);
            instanceData->days[i]->break_length = minute * 60;
            assert(success == 1);

            for (k = 0; k < 2; k++)
            {
                l++;
                if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
                {
                    fprintf(stderr,"Error skipping line %d.\n",l);
                    return SCIP_READERROR;
                }
            }
        }
        l++;
        if( !fgets(strbuffer, SCIP_MAXSTRLEN, inFILE) )
        {
            fprintf(stderr,"Error skipping line %d.\n",l);
            return SCIP_READERROR;
        }

        assert(!fgets(strbuffer, SCIP_MAXSTRLEN, inFILE));
        fclose(inFILE);
    }

   /* NB: ---preprocessing--- adjust the start date to the next work-day and the end date to the preceeding work-day, if
      necessary: */
   /* startDate is either a sunday, or a saturday but saturdays are work-free => move startDate to next monday */
   while( (((instanceData->startDate)->weekday == 5) && (! instanceData->workOnSaturdays)) || ((instanceData->startDate)->weekday == 6) )
   {
      SCIP_CALL( incrementDate(instanceData->startDate,1) );
   }
   /* endDate is either a sunday, or a saturday but saturdays are work-free => move endDate to the previous saturday or,
      if saturdays are work-free, friday */
   while( (((instanceData->endDate)->weekday == 5) && (! instanceData->workOnSaturdays)) || ((instanceData->endDate)->weekday == 6) )
   {
      SCIP_CALL( decrementDate(instanceData->endDate,1) );
   }


   return SCIP_OKAY;
}



/* set weekday for given date (by comparing with reference date macros) */
SCIP_RETCODE setWeekday(
   int day,                               /**< given day (1-31) */
   int month,                             /**< given month (1-12) */
   int year,                              /**< given year */
   int* weekday                           /**< pointer to integer in which to store the corresponding weekday index (0--6) */
   )
{
   int days_since_refdate  = 0;
   int days_per_month[12]  = {31,28,31,30,31,30,31,31,30,31,30,31};
   int ref_weekday = REFDATE_WEEKDAY;
   int i;

   for (i = REFDATE_YEAR; i < year; i++)
   {
      days_since_refdate += 365 + (int)((i % 4) == 0);
   }
   for (i = 0; i < month-1; i++)
   {
      days_since_refdate += days_per_month[i];
      if( (i == 1) && ((year % 4) == 0) ) /* february in a leap year */
         days_since_refdate++;
   }
   days_since_refdate += day-1;
   *weekday = (ref_weekday + days_since_refdate) % 7;

   return SCIP_OKAY;
}



/** initialize instance data struct */
SCIP_RETCODE initInstanceData(
   SCIP* scip,                            /**< SCIP pointer */
   instance_data* instanceData,           /**< pointer to struct in which instance data will be stored */
   int nC                                 /**< number of customers in data set (or -1 to init. everything with NULL) */
   )
{
   int k;

   if( nC <= 0 ) /* use -1 for any non-positive nC-value */
   {
      assert( instanceData != NULL );
      instanceData->nC = -1;
      instanceData->customerIDs            = NULL;
      instanceData->coordinates            = NULL;
      /* instanceData->timewindows         = NULL */
      instanceData->deliveryFreq           = NULL;
      instanceData->availability           = NULL;
      instanceData->estim_availability     = NULL;
      instanceData->appointment_sens       = NULL;
      instanceData->appointment_tol        = NULL;
      instanceData->t_service              = NULL;
      instanceData->estim_t_service        = NULL;
      instanceData->driver_bond            = NULL;
      instanceData->t_walk                 = NULL;
      instanceData->grade                  = NULL;
      instanceData->revenue                = NULL;
      instanceData->riskfactor             = NULL;
      /* instanceData->dependencies        = NULL; */
      instanceData->depot_coordinates      = NULL;
      instanceData->depot_out_distAM       = NULL;
      instanceData->depot_out_distNoon     = NULL;
      instanceData->depot_out_distPM       = NULL;
      instanceData->depot_out_t_travelAM   = NULL;
      instanceData->depot_out_t_travelNoon = NULL;
      instanceData->depot_out_t_travelPM   = NULL;
      instanceData->depot_in_distAM        = NULL;
      instanceData->depot_in_distNoon      = NULL;
      instanceData->depot_in_distPM        = NULL;
      instanceData->depot_in_t_travelAM    = NULL;
      instanceData->depot_in_t_travelNoon  = NULL;
      instanceData->depot_in_t_travelPM    = NULL;
      instanceData->workOnSaturdays        = FALSE;
      instanceData->startDate              = NULL;
      instanceData->endDate                = NULL;
   }
   else
   {
      assert( instanceData != NULL );
      instanceData->nC = nC;
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->customerIDs), instanceData->nC) );
      assert( instanceData->customerIDs != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->coordinates), instanceData->nC) );
      assert( instanceData->coordinates != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->timeWindows), instanceData->nC) );
      assert( instanceData->timeWindows != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->deliveryFreq), instanceData->nC) );
      assert( instanceData->deliveryFreq != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->availability), instanceData->nC) );
      assert( instanceData->availability != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->estim_availability), instanceData->nC) );
      assert( instanceData->estim_availability != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->appointment_sens), instanceData->nC) );
      assert( instanceData->appointment_sens != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->appointment_tol), instanceData->nC) );
      assert( instanceData->appointment_tol != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_service), instanceData->nC) );
      assert( instanceData->t_service != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->estim_t_service), instanceData->nC) );
      assert( instanceData->estim_t_service != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->driver_bond), instanceData->nC) );
      assert( instanceData->driver_bond != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_walk), instanceData->nC) );
      assert( instanceData->t_walk != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->grade), instanceData->nC) );
      assert( instanceData->grade != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->revenue), instanceData->nC) );
      assert( instanceData->revenue != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->riskfactor), instanceData->nC) );
      assert( instanceData->riskfactor != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distAM), instanceData->nC) );
      assert( instanceData->distAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distNoon), instanceData->nC) );
      assert( instanceData->distNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distPM), instanceData->nC) );
      assert( instanceData->distPM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelAM), instanceData->nC) );
      assert( instanceData->t_travelAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelNoon), instanceData->nC) );
      assert( instanceData->t_travelNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelPM), instanceData->nC) );
      assert( instanceData->t_travelPM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->lastServiceDate), instanceData->nC) );
      assert( instanceData->lastServiceDate != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_coordinates), 2) );
      assert( instanceData->depot_coordinates != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_distAM), instanceData->nC) );
      assert( instanceData->depot_out_distAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_distNoon), instanceData->nC) );
      assert( instanceData->depot_out_distNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_distPM), instanceData->nC) );
      assert( instanceData->depot_out_distPM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_t_travelAM), instanceData->nC) );
      assert( instanceData->depot_out_t_travelAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_t_travelNoon), instanceData->nC) );
      assert( instanceData->depot_out_t_travelNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_out_t_travelPM), instanceData->nC) );
      assert( instanceData->depot_out_t_travelPM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_distAM), instanceData->nC) );
      assert( instanceData->depot_in_distAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_distNoon), instanceData->nC) );
      assert( instanceData->depot_in_distNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_distPM), instanceData->nC) );
      assert( instanceData->depot_in_distPM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_t_travelAM), instanceData->nC) );
      assert( instanceData->depot_in_t_travelAM != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_t_travelNoon), instanceData->nC) );
      assert( instanceData->depot_in_t_travelNoon != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->depot_in_t_travelPM), instanceData->nC) );
      assert( instanceData->depot_in_t_travelPM != NULL );
      SCIP_CALL( SCIPallocBlockMemory(scip, &(instanceData->startDate)) );
      assert( instanceData->startDate != NULL );
      SCIP_CALL( SCIPallocBlockMemory(scip, &(instanceData->endDate)) );
      assert( instanceData->endDate != NULL );
      for (k = 0; k < instanceData->nC; k++)
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->coordinates[k]), 2) );
         assert( instanceData->coordinates[k] != NULL );
         instanceData->timeWindows[k] = NULL; /* list elements will be allocated on an as-needed basis */
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->estim_t_service[k]), 2) );
         assert( instanceData->estim_t_service[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distAM[k]), instanceData->nC) );
         assert( instanceData->distAM[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distNoon[k]), instanceData->nC) );
         assert( instanceData->distNoon[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->distPM[k]), instanceData->nC) );
         assert( instanceData->distPM[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelAM[k]), instanceData->nC) );
         assert( instanceData->t_travelAM[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelNoon[k]), instanceData->nC) );
         assert( instanceData->t_travelNoon[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(instanceData->t_travelPM[k]), instanceData->nC) );
         assert( instanceData->t_travelPM[k] != NULL );
         SCIP_CALL( SCIPallocBlockMemory(scip, &(instanceData->lastServiceDate[k])) );
         assert( instanceData->lastServiceDate[k] != NULL );
      }
      /* instanceData->timewindows ... */
      /* instanceData->dependencies ... */

   }

   return SCIP_OKAY;
}



/** de-initialize instance data struct (clean up memory) */
SCIP_RETCODE deinitInstanceData(
   SCIP* scip,                            /**< SCIP pointer */
   instance_data* instanceData            /**< pointer to struct in which instance data is stored */
   )
{
   int k;
   timeWindow* tmp_tw = NULL;

   for (k = (instanceData->nC)-1; k >= 0; k--)
   {
      SCIPfreeBlockMemory(scip, &(instanceData->lastServiceDate[k]));
      SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelPM[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelNoon[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelAM[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->distPM[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->distNoon[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->distAM[k]), instanceData->nC);
      SCIPfreeBlockMemoryArray(scip, &(instanceData->estim_t_service[k]), 2);
      while( instanceData->timeWindows[k] != NULL )
      {
         tmp_tw = instanceData->timeWindows[k];
         instanceData->timeWindows[k] = (instanceData->timeWindows[k])->next;
         SCIPfreeBlockMemory(scip, &(tmp_tw->twDate));
         SCIPfreeBlockMemory(scip, &tmp_tw);
      }
      SCIPfreeBlockMemoryArray(scip, &(instanceData->coordinates[k]), 2);
   }
   SCIPfreeBlockMemory(scip, &(instanceData->endDate));
   SCIPfreeBlockMemory(scip, &(instanceData->startDate));
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_t_travelPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_t_travelNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_t_travelAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_distPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_distNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_in_distAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_t_travelPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_t_travelNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_t_travelAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_distPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_distNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_out_distAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->depot_coordinates), 2);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->lastServiceDate), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->t_travelAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->distPM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->distNoon), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->distAM), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->riskfactor), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->revenue), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->grade), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->t_walk), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->driver_bond), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->estim_t_service), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->t_service), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->appointment_tol), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->appointment_sens), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->estim_availability), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->availability), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->deliveryFreq), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->timeWindows), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->coordinates), instanceData->nC);
   SCIPfreeBlockMemoryArray(scip, &(instanceData->customerIDs), instanceData->nC);
   if(instanceData->days != NULL)
   {
       SCIPfreeBlockMemoryArray(scip, &(instanceData->days), instanceData->nDays);
   }
   /* instanceData->dependencies ... */

   return SCIP_OKAY;
}



/** create (initialize and fill) model data struct */
SCIP_RETCODE createModelData(
   SCIP* scip,                            /**< SCIP pointer */
   model_data* modelData,                 /**< pointer to struct in which model data will be stored */
   instance_data* instanceData            /**< pointer to struct in which instance data is stored */
   )
{
   int i,j,k;

   assert( modelData != NULL );
   assert( instanceData != NULL );
   /* initialization and memory allocation */
   modelData->nC               = instanceData->nC + 1;   /* number of customers + 1 for the depot */
   modelData->shift_start      = DEFAULT_SHIFT_START;    /* use default 06:00 for now; could later be adapted to possible additional data read from input */
   modelData->shift_end        = DEFAULT_SHIFT_END;      /* use default 22:00 for now; -"- */
   modelData->maxDelayEvents   = DEFAULT_MAXDELAYEVENTS; /* use default 10 for now; -"- */
   modelData->nDays            = instanceData->nDays;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->customerIDs), modelData->nC) );
   assert( modelData->customerIDs != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->obj), modelData->nC) );
   assert( modelData->obj != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->timeWindows), modelData->nC) );
   assert( modelData->timeWindows != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_service), modelData->nC) );
   assert( modelData->t_service != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_service_maxDev), modelData->nC) );
   assert( modelData->t_service_maxDev != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelAM), modelData->nC) );
   assert( modelData->t_travelAM != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelNoon), modelData->nC) );
   assert( modelData->t_travelNoon != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelPM), modelData->nC) );
   assert( modelData->t_travelPM != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travel), modelData->nC) );
   assert( modelData->t_travel != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travel_maxDev), modelData->nC) );
   assert( modelData->t_travel_maxDev != NULL );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(modelData->startDate)) );
   assert( modelData->startDate != NULL );
   SCIP_CALL( SCIPallocBlockMemory(scip, &(modelData->endDate)) );
   assert( modelData->endDate != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->days), modelData->nDays) );
   assert( modelData->days != NULL );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->day_sizes), modelData->nDays) );
   assert( modelData->days != NULL );
   for (k = 0; k < (modelData->nC); k++)
   {
      modelData->timeWindows[k] = NULL; /* list elements will be allocated on an as-needed basis */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelAM[k]), modelData->nC) );
      assert( modelData->t_travelAM[k] != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelNoon[k]), modelData->nC) );
      assert( modelData->t_travelNoon[k] != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travelPM[k]), modelData->nC) );
      assert( modelData->t_travelPM[k] != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travel[k]), modelData->nC) );
      assert( modelData->t_travel[k] != NULL );
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(modelData->t_travel_maxDev[k]), modelData->nC) );
      assert( modelData->t_travel_maxDev[k] != NULL );
   }

   /* derive and/or copy data parameters from instance data */
   modelData->workOnSaturdays      = instanceData->workOnSaturdays; /* kept for information purposes: are saturdays included as workdays? (TRUE/FALSE) */
   (modelData->startDate)->day     = (instanceData->startDate)->day;
   (modelData->startDate)->month   = (instanceData->startDate)->month;
   (modelData->startDate)->year    = (instanceData->startDate)->year;
   (modelData->startDate)->weekday = (instanceData->startDate)->weekday;
   (modelData->endDate)->day       = (instanceData->endDate)->day;
   (modelData->endDate)->month     = (instanceData->endDate)->month;
   (modelData->endDate)->year      = (instanceData->endDate)->year;
   (modelData->endDate)->weekday   = (instanceData->endDate)->weekday;
   for (i = 0; i < (modelData->nC); i++)
   {
      /* copy/compute travel time matrices: */
      for (j = 0; j < (modelData->nC); j++)
      {
         if( i < (modelData->nC)-1 ) /* "from" customer is actual customer */
         {
            if( j < (modelData->nC)-1 ) /* "to" customer is actual customer */
            {
               modelData->t_travelAM[i][j]   = instanceData->t_travelAM[i][j];
               modelData->t_travelNoon[i][j] = instanceData->t_travelNoon[i][j];
               modelData->t_travelPM[i][j]   = instanceData->t_travelPM[i][j];
               modelData->t_travel[i][j]     = (int)ceil((double)(instanceData->t_travelAM[i][j]+instanceData->t_travelNoon[i][j]+instanceData->t_travelPM[i][j])/3.0);
            }
            else /* "to" customer is the depot */
            {
               modelData->t_travelAM[i][j]   = instanceData->depot_in_t_travelAM[i];
               modelData->t_travelNoon[i][j] = instanceData->depot_in_t_travelNoon[i];
               modelData->t_travelPM[i][j]   = instanceData->depot_in_t_travelPM[i];
               modelData->t_travel[i][j]     = (int)ceil((double)(instanceData->depot_in_t_travelAM[i]+instanceData->depot_in_t_travelNoon[i]+instanceData->depot_in_t_travelPM[i])/3.0);
            }
         }
         else /* "from" customer is the depot */
         {
            if( j < (modelData->nC)-1 ) /* "to" customer is actual customer */
            {
               modelData->t_travelAM[i][j]   = instanceData->depot_out_t_travelAM[j];
               modelData->t_travelNoon[i][j] = instanceData->depot_out_t_travelNoon[j];
               modelData->t_travelPM[i][j]   = instanceData->depot_out_t_travelPM[j];
               modelData->t_travel[i][j]     = (int)ceil((double)(instanceData->depot_out_t_travelAM[j]+instanceData->depot_out_t_travelNoon[j]+instanceData->depot_out_t_travelPM[j])/3.0);
            }
            else /* "to" customer is the depot */
            {
               assert( (i == modelData->nC-1) && (j == modelData->nC-1) );
               modelData->t_travelAM[i][j]   = 0.0;
               modelData->t_travelNoon[i][j] = 0.0;
               modelData->t_travelPM[i][j]   = 0.0;
               modelData->t_travel[i][j]     = 0.0;
            }
         }
         k = (int)MAX3( modelData->t_travelAM[i][j], modelData->t_travelNoon[i][j], modelData->t_travelPM[i][j] ); /* MAX3 macro (maximum of 3 values) is defined by SCIP */
         modelData->t_travel_maxDev[i][j] = k - modelData->t_travel[i][j];
      }
      /* copy/compute customerIDs, service times and assoc. max. deviations, objective function coefficients and lists of delivery time windows */
      if( i < (modelData->nC)-1 )
      {
         modelData->customerIDs[i] = instanceData->customerIDs[i];
         /* determine (average) service times and assoc. max. deviations; for details, see function setModelServiceTimes() */
         SCIP_CALL( setModelServiceTimes(&(modelData->t_service[i]), &(modelData->t_service_maxDev[i]), instanceData->t_walk[i], instanceData->t_service[i], instanceData->estim_t_service[i]) );
         /* determine objective function coefficients; for details, see function setModelObjectiveCoeff() */
         SCIP_CALL( setModelObjectiveCoeff(&(modelData->obj[i]),i,instanceData) );
         /* determine delivery time windows (taking concrete-date deliveries as well as recurring deliveries
            (wr.t. frequency and last delivery date) into account w.r.t. whole planning period);
            for details, see function setModelTimeWindows() */
         SCIP_CALL( setModelTimeWindows(scip,&(modelData->timeWindows[i]),i,instanceData) );
      }
      else
      {
         modelData->customerIDs[i]      = -1;   /* depot */
         modelData->t_service[i]        = 0;
         modelData->t_service_maxDev[i] = 0;
         modelData->timeWindows[i]      = NULL; /* leaving this empty for now -- start and end time of shifts are stored separately anyway; could become useful in possible later extensions */
         modelData->obj[i]              = 0.0;
      }
   }
   /* day properties */
   if(instanceData->days == NULL) /* set default data */
   {
       for(i = 0; i < modelData->nDays; i++)
       {
           SCIP_CALL( SCIPallocBlockMemory(scip, &(modelData->days[i])) );
           modelData->days[i]->index = i;
           modelData->days[i]->shift_start =    DEFAULT_SHIFT_START;
           modelData->days[i]->shift_end =      DEFAULT_SHIFT_END;
           modelData->days[i]->break_start =    DEFAULT_BREAK_START;
           modelData->days[i]->break_end =      DEFAULT_BREAK_END;
           modelData->days[i]->break_length =   DEFAULT_BREAK_LENGTH;
       }
   }
   else /* copy from instance data */
   {
       for(i = 0; i < modelData->nDays; i++)
       {
           SCIP_CALL( SCIPduplicateMemory(scip, &(modelData->days[i]), instanceData->days[i]) );
       }
   }
   modelData->neighbors = NULL; /* when being created, neighbor information is not yet present (corresp. to full digraph for every day); will be filled in preprocessing eventually */
   /* default day sizes */
   for(i = 0; i < modelData->nDays; i++)
   {
       modelData->day_sizes[i] = modelData->nC - 1;
   }

   return SCIP_OKAY;
}



/** de-initialize model data struct (clean up memory) */
SCIP_RETCODE deinitModelData(
   SCIP* scip,                            /**< SCIP pointer */
   model_data* modelData                  /**< pointer to struct in which model data is stored */
   )
{
   int j,k;
   modelWindow* tmp_tw = NULL;
   neighbor* tmp_nb    = NULL;

   for (k = (modelData->nC)-1; k >= 0; k--)
   {
      if( modelData->neighbors != NULL )
      {
          if(modelData->neighbors[k] != NULL) {
              for (j = modelData->nDays - 1; j >= 0; j--) {
                  if (modelData->neighbors[k][j] != NULL) {
                      while (modelData->neighbors[k][j] != NULL) {
                          tmp_nb = modelData->neighbors[k][j];
                          modelData->neighbors[k][j] = (modelData->neighbors[k][j])->next;
                          SCIPfreeBlockMemory(scip, &tmp_nb);
                      }
                  }
              }
              if(modelData->neighbors[k] != NULL)
              {
                  SCIPfreeBlockMemoryArray(scip, &(modelData->neighbors[k]), modelData->nDays);
              }
          }
      }
      SCIPfreeBlockMemoryArray(scip, &(modelData->t_travel_maxDev[k]), modelData->nC);
      SCIPfreeBlockMemoryArray(scip, &(modelData->t_travel[k]), modelData->nC);
      SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelPM[k]), modelData->nC);
      SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelNoon[k]), modelData->nC);
      SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelAM[k]), modelData->nC);
      while( modelData->timeWindows[k] != NULL )
      {
         tmp_tw = modelData->timeWindows[k];
         modelData->timeWindows[k] = (modelData->timeWindows[k])->next;
         SCIPfreeBlockMemory(scip, &tmp_tw);
      }
   }
   if( modelData->neighbors != NULL )
   {
      SCIPfreeBlockMemoryArray(scip, &(modelData->neighbors), modelData->nC);
   }
   SCIPfreeBlockMemoryArray(scip, &(modelData->day_sizes), modelData->nDays);
   for(k = 0; k < modelData->nDays; k++)
   {
       if(modelData->days[k] != NULL)
       {
           SCIPfreeBlockMemory(scip, &(modelData->days[k]));
       }
   }
   SCIPfreeBlockMemoryArray(scip, &(modelData->days), modelData->nDays);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_service_maxDev), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_service), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_travel_maxDev), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_travel), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelPM), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelNoon), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->t_travelAM), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->timeWindows), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->obj), modelData->nC);
   SCIPfreeBlockMemoryArray(scip, &(modelData->customerIDs), modelData->nC);
   SCIPfreeBlockMemory(scip, &(modelData->startDate));
   SCIPfreeBlockMemory(scip, &(modelData->endDate));
   return SCIP_OKAY;
}



/** derive and set model service times and assoc. max. deviations from instance data */
SCIP_RETCODE setModelServiceTimes(
   int* t_service,                        /**< pointer to where model service time is to be stored */
   int* t_service_maxDev,                 /**< pointer to where model max. service time deviation is to be stored */
   int twalk,                             /**< instance data walking time at customer */
   int ts,                                /**< instance data service time */
   int* ets                               /**< instance data estimated service time interval (2-array) */
   )
{
   /* Note: This function can be modified without having to change anything else, if the service times and
      max. deviation values are to be computed differently than in this first suggestion. */
   /* instance data provides a service time ts and an estimated service time interval [ets[0],ets[1]] (not necessarily
      containing ts) along with the time for walking (from truck to customer). In case ts lies in the interval, we take
      ts as the model service time and the difference to ets[1] as the max. deviation; otherwise, model service time
      becomes 0.5*(min(max(ts,1),ets[0])+max(ts,ets[1])) and the max. deviation becomes the difference from model
      service time to max(ts,ets[1]). Moreover, walking time will be integrated into (added to) service times. */
   /* NB: instance data contains service and walking times in minutes;
      these are converted to seconds for model data here. */
   assert( t_service != NULL );
   assert( t_service_maxDev != NULL );
   assert( ets != NULL );
   *t_service        = 60*ets[0];
   *t_service_maxDev = 60*(ets[1]-ets[0]);

   *t_service        += 60*twalk;
   return SCIP_OKAY;
}



/** derive and set model objective function coefficient for customer i from instance data */
SCIP_RETCODE setModelObjectiveCoeff(
   SCIP_Real* obj,                        /**< pointer to where i-th model objective coefficient is to be stored */
   int  i,                                /**< customer index for which to derive objective coefficient */
   instance_data* instanceData            /**< instance data struct */
   )
{
   SCIP_Real estim_availability, appointment_sens, grade, revenue, riskfactor;
   /* Note: This function can be modified without having to change anything else, if the objective coefficient values
      are to be computed differently than in this first suggestion. */
   /* instance data contains several parameters that are taken into account to derive objective function coefficients
      that ultimately determine how much lateness at a customer is being penalized: availability, estim_availability,
      appointment_sens, appointment_tol, driver_bond, grade, revenue, riskfactor */
   assert( obj != NULL );
   assert( instanceData != NULL );
   assert( (i >= 0) && (i < instanceData->nC) );

   /* for now, we neglect "availabliity" (binary values), "driver_bond" ([0,1]), and "appointment_tol" (rescheduling
      tolerance, {0,1,..,20}), since it is entirely unclear how these data would influence a penalty on service
      lateness */

   /* normalize values to numbers in [0,1] and such that higher value leads to higher penalty */
   estim_availability = instanceData->estim_availability[i];                    /* [0,1] estim. probabilities */
   appointment_sens   = 1.0/(SCIP_Real)(instanceData->appointment_sens[i]);     /* 1,2,3 (A,B,C) appointment sensitivity */
   grade              = 1.0/(SCIP_Real)(instanceData->grade[i]);                /* 1--8 (A--D3) customer grade/value */
   revenue            = (SCIP_Real)(instanceData->revenue[i])/5000.0;           /* [0,5000] customer revenue values in EUR */
   riskfactor         = instanceData->riskfactor[i];                            /* [0,1] drop out (if unsatisfied) risk/probability */
   /* for now, assume order of importance: grade, revenue, availability, appointment sensitivity, drop-out risk */
   *obj = 1000*grade + 100*revenue + 10*estim_availability + 1*appointment_sens + 0.1*riskfactor;

   return SCIP_OKAY;
}



/** derive and set model delivery time windows for customer i from instance data */
SCIP_RETCODE setModelTimeWindows(
   SCIP* scip,                            /**< SCIP pointer */
   modelWindow** tw,                      /**< pointer to where i-th list of model delivery time windows is to be stored */
   int  i,                                /**< customer index for which to derive time windows */
   instance_data* instanceData            /**< instance data struct */
   )
{
   modelWindow* tmp_tw = NULL;
   modelWindow* twHead = NULL;
   timeWindow*  inst_tw = NULL;
   int deliveryFreq;
   date* lastServiceDate;
   date* inst_tw_date = NULL;
   SCIP_Bool inst_tw_recurring;
   int start_t, end_t;
   int itmp = 0;
   int itmp2 = 0;
   int weeklyWorkdays = 0;
   int endDateDay = 0;

   assert( *tw == NULL ); /* list should have been initialized with NULL */
   assert( instanceData != NULL );
   assert( (i >= 0) && (i < instanceData->nC) );
   inst_tw = instanceData->timeWindows[i];
   if( inst_tw == NULL )
   {
      fprintf(stderr,"time window list of customer %d is empty.\n",i);
      return SCIP_ERROR;
   }
   assert( inst_tw != NULL );

   weeklyWorkdays = ( instanceData->workOnSaturdays ) ? 6 : 5;
   SCIP_CALL( setModelDay(&endDateDay,instanceData->endDate,instanceData) ); /* model day-index of end date of planning period */
   deliveryFreq    = instanceData->deliveryFreq[i];
   lastServiceDate = instanceData->lastServiceDate[i];
   assert( lastServiceDate != NULL );
   /* In case no last service date was given but the given delivery frequency is different from 1 (which would allow to
      schedule delivery in every week of the planning period), the deliveryFreq is reset to 1 (because other cases
      including 0 do not really make sense): */
   if( (lastServiceDate->weekday == -1) && (deliveryFreq != 1) )
   {
#ifdef SCIP_DEBUG
      SCIPinfoMessage(scip,stdout,"[Warning:] delivery frequency of customer %3d reset to 1 (was %d) because no associated last service date was provided.\n",i,deliveryFreq);
#endif
      instanceData->deliveryFreq[i] = 1;
      deliveryFreq = 1;
   }

   twHead = *tw;
   assert( twHead == NULL );
   /* traverse instance data time-window list and create time-window list in model data: */
   while( inst_tw != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
      inst_tw_recurring = inst_tw->recurring;
      inst_tw_date      = inst_tw->twDate;
      assert( inst_tw_date != NULL );
      start_t = 60*(inst_tw->start_min + 60*inst_tw->start_hour);
      end_t   = 60*(inst_tw->end_min + 60*inst_tw->end_hour);
      if( inst_tw_recurring ) /* recurring time windows => every resp. weekday over whole planning period, subject to deliveryFreq w.r.t. lastServiceDate */
      {

         assert( tmp_tw != NULL );
         /* set tmp_tw->day to model day-index of first specified weekday after start of planning period: */
         SCIP_CALL( setModelDay(&(itmp),inst_tw_date,instanceData) );
         if(itmp > endDateDay)
         {
             inst_tw = inst_tw->next;
             SCIPfreeBlockMemory(scip, &tmp_tw);
             continue;
         }
         assert( itmp <= endDateDay );
         if( lastServiceDate->weekday == -1 || deliveryFreq == 0) /* no last service date given */
         {
            while( itmp <= endDateDay )
            {

               tmp_tw->day = itmp;
               tmp_tw->start_t = start_t;
               tmp_tw->end_t   = end_t;
               tmp_tw->next    = NULL;
               tmp_tw->weigth = inst_tw->weigth;
               if( twHead == NULL ) /* first list element */
               {
                  twHead = tmp_tw;
                  *tw = tmp_tw;
                  tmp_tw = tmp_tw->next;
                  assert( tw != NULL );
               }
               else /* append new list element to end of the current list */
               {
                  twHead->next = tmp_tw;
                  twHead       = twHead->next;
                  tmp_tw       = tmp_tw->next;
               }
               assert( twHead != NULL );
               assert( tmp_tw == NULL );
               itmp += weeklyWorkdays;
               if( itmp <= endDateDay )
               {
                  SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
                  assert( tmp_tw != NULL );
               }
            }
         }
         else if( deliveryFreq == 0 ) /* a last service date was given, but deliveryFreq calls for one-time delivery */
         {
            /* in this case, we ignore the last serice date (for now) and proceed as if none was given, setting the time
               window to the first available one on the specified weekday */
            assert( lastServiceDate->weekday != -1 );
            tmp_tw->day = itmp;
            tmp_tw->start_t = start_t;
            tmp_tw->end_t   = end_t;
            tmp_tw->next    = NULL;
            tmp_tw->weigth = inst_tw->weigth;
            if( twHead == NULL ) /* first list element */
            {
               twHead = tmp_tw;
               *tw = tmp_tw;
               tmp_tw = tmp_tw->next;
               assert( tw != NULL );
            }
            else /* append new list element to end of the current list */
            {
               twHead->next = tmp_tw;
               twHead       = twHead->next;
               tmp_tw       = tmp_tw->next;
            }
            assert( twHead != NULL );
            assert( tmp_tw == NULL );
         }
         else /* a last service date was given and delivery shall only be allowed in a deliveryFreq-weekly rhythm
                 (starting from last service date) over the planning period; we create the time windows for the weeks in
                 deliveryFreq week-intervals from the week in which the last service occurred. */
            /* TO DO: the model is currently not adapted to having multiple deliveries to the same customer over the
               planning period- This could be handled by creating customer copies, or via additional constraints and
               adjusted "visit each customer once"-constraints etc. (Alternatively, for now, all but the first time
               window correponding to the weekday under consideration could be ignored when setting up the model.) */
         {
            assert( deliveryFreq > 0 );
            assert( lastServiceDate->weekday != -1 );
            assert( (lastServiceDate->day != -1) && (lastServiceDate->month != -1) && (lastServiceDate->year != -1) );
            itmp2 = -getWorkdayDiff(lastServiceDate,instanceData->startDate,instanceData); /* model day of lastServiceDate is itmp2<0 */
            /* must ensure itmp, the model-day of first given weekday after planning period start date, is in a week
               that is a multiple of deliveryFreq weeks after the week containing the lastServiceDate
               (it can happen that last service date is longer ago than the deliveryFreq would suggest -> in this case
               we must move itmp up to the first week adhering to deliveryFreq after the startDate) */
            while( itmp2 + deliveryFreq*weeklyWorkdays +(weeklyWorkdays-1-lastServiceDate->weekday) < itmp )
            {
               /* the whole week where the last service date lies is longer ago (w.r.t. day itmp under consideration)
                  than deliveryFreq weeks => increment itmp2 by deliveryFreq*weeklyWorkdays */
               itmp2 += deliveryFreq*weeklyWorkdays;
            }
            if( (itmp2 + deliveryFreq*weeklyWorkdays -lastServiceDate->weekday <= itmp) && (itmp2 + deliveryFreq*weeklyWorkdays +(weeklyWorkdays-1-lastServiceDate->weekday) >= itmp) )
            {
               /* itmp falls into a week that is a multiple of deliveryFreq weeks from the one where last service
                  took place -> start here (first time window associated with weekday under consideration) */
               while( itmp <= endDateDay )
               {
                  tmp_tw->day = itmp;
                  tmp_tw->start_t = start_t;
                  tmp_tw->end_t   = end_t;
                  tmp_tw->next    = NULL;
                  tmp_tw->weigth = inst_tw->weigth;
                  if( twHead == NULL ) /* first list element */
                  {
                     twHead = tmp_tw;
                     *tw = tmp_tw;
                     tmp_tw = tmp_tw->next;
                     assert( tw != NULL );
                  }
                  else /* append new list element to end of the current list */
                  {
                     twHead->next = tmp_tw;
                     twHead       = twHead->next;
                     tmp_tw       = tmp_tw->next;
                  }
                  assert( twHead != NULL );
                  assert( tmp_tw == NULL );
                  itmp += deliveryFreq*weeklyWorkdays;
                  if( itmp <= endDateDay )
                  {
                     SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
                     assert( tmp_tw != NULL );
                  }
               }
            }
            else /* itmp is too early -> move first time window to week corresponding to deliveryFreq weeks after last service date */
            {
               assert( itmp2 + deliveryFreq*weeklyWorkdays -lastServiceDate->weekday > itmp );
               while( itmp2 + deliveryFreq*weeklyWorkdays -lastServiceDate->weekday > itmp )
               {
                  itmp += weeklyWorkdays;
                  if( itmp > endDateDay )
                  {
                     #ifdef SCIP_DEBUG
                     SCIPinfoMessage(scip,stdout,"[Warning:] a time window was moved outside of planning period due to associated delivery frequency and last service date -> not considered.\n");
                     #endif
                     /* printf("           details: lastServiceDate: %d.%d.%d (%d), deliveryFreq: %d, startDate: %d.%d.%d (%d), weekday of time window: %d, endDate: %d.%d.%d (%d)\n",
                        lastServiceDate->day,lastServiceDate->month,lastServiceDate->year,lastServiceDate->weekday,deliveryFreq,instanceData->startDate->day,instanceData->startDate->month,instanceData->startDate->year,instanceData->startDate->weekday,inst_tw_date->weekday,instanceData->endDate->day,instanceData->endDate->month,instanceData->endDate->year,instanceData->endDate->weekday); */
                     SCIPfreeBlockMemory(scip, &tmp_tw);
                     assert( tmp_tw == NULL );
                     if( itmp2 + deliveryFreq*weeklyWorkdays-lastServiceDate->weekday > itmp )
                     {
                        itmp = itmp2 + deliveryFreq*weeklyWorkdays-lastServiceDate->weekday;
                     }
                  }
               }
               assert( itmp2 + deliveryFreq*weeklyWorkdays +(weeklyWorkdays-1-lastServiceDate->weekday) >= itmp );
               while( itmp <= endDateDay )
               {
                  tmp_tw->day = itmp;
                  tmp_tw->start_t = start_t;
                  tmp_tw->end_t   = end_t;
                  tmp_tw->next    = NULL;
                  tmp_tw->weigth = inst_tw->weigth;
                  if( twHead == NULL ) /* first list element */
                  {
                     twHead = tmp_tw;
                     *tw = tmp_tw;
                     tmp_tw = tmp_tw->next;
                     assert( tw != NULL );
                  }
                  else /* append new list element to end of the current list */
                  {
                     twHead->next = tmp_tw;
                     twHead       = twHead->next;
                     tmp_tw       = tmp_tw->next;
                  }
                  assert( twHead != NULL );
                  assert( tmp_tw == NULL );
                  itmp += deliveryFreq*weeklyWorkdays;
                  if( itmp <= endDateDay )
                  {
                     SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
                     assert( tmp_tw != NULL );
                  }
               }
            }

         }
      }
      else /* non-recurring time windows => specific date range only (deliveryFreq, lastServiceDate are unused/irrelevant) */
      {
         assert( tmp_tw != NULL );
         SCIP_CALL( setModelDay(&(tmp_tw->day),inst_tw_date,instanceData) );
          if(tmp_tw->day > endDateDay)
          {
              inst_tw = inst_tw->next;
              SCIPfreeBlockMemory(scip, &tmp_tw);
              continue;
          }
         tmp_tw->start_t = start_t;
         tmp_tw->end_t   = end_t;
         tmp_tw->next    = NULL;
         tmp_tw->weigth = inst_tw->weigth;
         if( twHead == NULL ) /* first list element */
         {
            twHead = tmp_tw;
            *tw = tmp_tw;
            tmp_tw = tmp_tw->next;
            assert( tmp_tw == NULL );
            assert( *tw != NULL );
         }
         else /* append new list element to end of the current list, possibly after creating and appending elements for the days between current and first day */
         {
            itmp        = tmp_tw->day; /* last day in range */
            itmp2       = twHead->day+1;
            tmp_tw->day = itmp2;
            while( tmp_tw->day < itmp ) /* day range spans more than two consecutive days -> create and append elements for middle days */
            {
               twHead->next = tmp_tw;
               twHead       = twHead->next;
               tmp_tw       = tmp_tw->next;
               assert( tmp_tw == NULL );
               assert( twHead != NULL );
               SCIP_CALL( SCIPallocBlockMemory(scip, &tmp_tw) );
               assert( tmp_tw != NULL );
               tmp_tw->day     = ++itmp2;
               tmp_tw->start_t = start_t;
               tmp_tw->end_t   = end_t;
               tmp_tw->next    = NULL;
               tmp_tw->weigth = inst_tw->weigth;
            }
            twHead->next = tmp_tw;
            twHead       = twHead->next;
            tmp_tw       = tmp_tw->next;
            assert( tmp_tw == NULL );
            assert( twHead != NULL );
         }
      }

      inst_tw = inst_tw->next;
   }
   return SCIP_OKAY;
}



/* increment the date in a date struct by a given number of days */
SCIP_RETCODE incrementDate(
   date* d,                               /**< date structure */
   int   k                                /**< number of days to increase date d by */
   )
{
   int days_per_month[12]  = {31,28,31,30,31,30,31,31,30,31,30,31};

   if( d->day < days_per_month[d->month-1] ) /* only day counter needs to be increased */
   {
      d->day += 1;
   }
   else /* date increase brings us into next month: */
   {
      assert( d->day == days_per_month[d->month-1] );
      d->day = 1;
      if( d->month < 12 ) /* just increase month counter */
      {
         d->month += 1;
      }
      else /* month increase brings us into next year: */
      {
         assert( d->month == 12 );
         d->month = 1;
         d->year += 1;
      }
   }
   d->weekday = (d->weekday + 1) % 7;

   return SCIP_OKAY;
}



/* decrement the date in a date struct by a given number of days */
SCIP_RETCODE decrementDate(
   date* d,                               /**< date structure */
   int   k                                /**< number of days to decrease date d by */
   )
{
   int days_per_month[12]  = {31,28,31,30,31,30,31,31,30,31,30,31};

   if( d->day > 1 ) /* only day counter needs to be decreased */
   {
      d->day -= 1;
   }
   else /* date decrease brings us into previous month: */
   {
      assert( d->day == 1 );
      if( d->month > 1 ) /* just decrease month counter */
      {
         d->month -= 1;
      }
      else /* month decrease brings us into previous year: */
      {
         d->month = 12;
         d->year -= 1;
      }
      d->day = days_per_month[d->month-1];
   }
   d->weekday = (d->weekday == 0) ? 6 : (d->weekday - 1);

   return SCIP_OKAY;
}



/* set model time window (work)day index */
SCIP_RETCODE setModelDay(
   int* day,                              /* pointer to where to store workday index */
   date* twDate,                          /* date struct of the day under consideration */
   instance_data* instanceData            /* instance data struct */
   )
{

   date* tmpDate = (date* )malloc(1*sizeof(date));
   assert( tmpDate != NULL );
   tmpDate->day     = (instanceData->startDate)->day;
   tmpDate->month   = (instanceData->startDate)->month;
   tmpDate->year    = (instanceData->startDate)->year;
   tmpDate->weekday = (instanceData->startDate)->weekday;

   if( twDate->day != -1 ) /* the case that the day under consideration is a concrete date */
   {
      (*day) = 0;
      while( (tmpDate->day != twDate->day) || (tmpDate->month != twDate->month) || (tmpDate->year != twDate->year) )
      {
         SCIP_CALL( incrementDate(tmpDate,1) );
         if( (tmpDate->weekday < 5) || ((tmpDate->weekday == 5) && (instanceData->workOnSaturdays)) )
         {
            (*day)++;
         }
      }
   }
   else /* the case that the day under consideration is given only as a weekday */
   {
      assert( (twDate->day==-1) && (twDate->month==-1) && (twDate->year==-1) );
      assert( twDate->weekday >= 0 );
      (instanceData->workOnSaturdays) ? assert( twDate->weekday <= 5 ) : assert( twDate->weekday < 5 );
      (*day) = 0;
      if( twDate->weekday >= tmpDate->weekday ) /* first day with given weekday lies in same week as start of planning period */
      {
         /* start of planning period = day 0 in the model => given weekday first occurs at this day: */
         (*day) = twDate->weekday - tmpDate->weekday;
         assert( *day >= 0 );
      }
      else /* first day with given weekday lies in the week after that where planning period starts */
      {
         /* given weekday is "model day -(tmpDate->weekday-twDate->weekday)" (<0)
            => +5 or 6 workdays gives the corresponding model day in next week */
         (*day) = (instanceData->workOnSaturdays) ? twDate->weekday - tmpDate->weekday + 6 : twDate->weekday - tmpDate->weekday + 5;
      }
   }

   free(tmpDate);

   return SCIP_OKAY;
}



/* get workday difference for given date pair */
int getWorkdayDiff(
   date* date1,                           /* date struct of the EARLIER date */
   date* date2,                           /* date struct of the LATER date */
   instance_data* instanceData            /* instance data struct */
   )
{
   int days;
   date* tmpDate = (date* )malloc(1*sizeof(date));
   assert( tmpDate != NULL );
   assert( (date1->day != -1) && (date1->month != -1) && (date1->year != -1) && (date1->weekday != -1) );
   assert( (date2->day != -1) && (date2->month != -1) && (date2->year != -1) && (date2->weekday != -1) );
   tmpDate->day     = (date1)->day;
   tmpDate->month   = (date1)->month;
   tmpDate->year    = (date1)->year;
   tmpDate->weekday = (date1)->weekday;

   days = 0;
   while( (tmpDate->day != date2->day) || (tmpDate->month != date2->month) || (tmpDate->year != date2->year) )
   {
      SCIP_CALL( incrementDate(tmpDate,1) );
      if( (tmpDate->weekday < 5) || ((tmpDate->weekday == 5) && (instanceData->workOnSaturdays)) )
      {
         days++;
      }
   }

   free(tmpDate);

   return days;
}



/* write model data to file */
SCIP_RETCODE writeModelData(
   model_data* modelData,                 /**< pointer to model data struct */
   const char* outputFile                 /**< output file name */
   )
{
   FILE* outFILE = NULL;
   int i,j,k;
   modelWindow* tmptw = NULL;
   neighbor* tmpnb    = NULL;

   assert( modelData != NULL );
   if( outputFile == NULL )
   {
      #ifdef SCIP_DEBUG
      printf("[Warning:] no output file name provided, aborting writeModelData.\n");
      #endif
      return SCIP_OKAY;
   }

   outFILE = fopen(outputFile,"w");
   assert( outFILE != NULL );

   /* int                   nC;                     /\**< number of customers + 1 (depot stored as last "customer") *\/ */
   fprintf(outFILE,"# number of customers + 1 (depot stored as last \"customer\") [int nC]\n");
   fprintf(outFILE,"%d\n",modelData->nC);

   /* int*                  customerIDs;            /\**< nC-array with customer IDs (not needed for model, but kept for output/information purposes; depot gets ID -1) *\/ */
   fprintf(outFILE,"# customer IDs (depot gets ID -1) [int* customerIDs]\n");
   assert( modelData->customerIDs != NULL );
   for( i = 0; i < modelData->nC; i++ )
   {
      fprintf(outFILE,"%d ",modelData->customerIDs[i]);
   }
   fprintf(outFILE,"\n");

   /* modelWindow**         timeWindows;            /\**< nC-array with linked list structures to store customers' delivery time windows *\/ */
   fprintf(outFILE,"# time window lists for each real customer (number of elements followed by elements' contents (day,start_t,end_t,weight)) [modelWindow** timeWindows]\n");
   for( i = 0; i < modelData->nC-1; i++ ) /* depot has NULL ... -> omitted */
   {
      tmptw = modelData->timeWindows[i];
      /* assert( tmptw != NULL ); */ /* USE ONLY IF assumption is that every customer actually has a time window... (it
         can happen (e.g. customer 32 in hessen data set) that last service date, delivery frequency and start/end of
         planning period yield an empty time window => customer is unreachable in planning period => can be discarded in
         preprocessing...). TO DO @logiball: this probably should not happen at all */
      k = 0;
      while( tmptw != NULL )
      {
         tmptw = tmptw->next;
         k++;
      }
      fprintf(outFILE,"%d ",k);
      tmptw = modelData->timeWindows[i];
      while( tmptw != NULL )
      {
         fprintf(outFILE,"%d %d %d %.2f ",tmptw->day,tmptw->start_t,tmptw->end_t, tmptw->weigth);
         tmptw = tmptw->next;
      }
      fprintf(outFILE,"\n");
   }

   /* SCIP_Real*            obj;                    /\**< nC-array with objective function coefficients for each customer *\/ */
   fprintf(outFILE,"# objective function coefficients [SCIP_Real* obj]\n");
   assert( modelData->obj != NULL );
   for( i = 0; i < modelData->nC; i++ )
   {
      fprintf(outFILE,"%f ",modelData->obj[i]);
   }
   fprintf(outFILE,"\n");

   /* int                   shift_start;            /\**< time when driver's daily work shift starts (no departure from depot before then) (in seconds) *\/ */
   fprintf(outFILE,"# shift start time [int shift_start]\n");
   fprintf(outFILE,"%d\n",modelData->shift_start);

   /* int                   shift_end;              /\**< time when driver's daily work shift ends (must return to depot by then) (in seconds) *\/ */
   fprintf(outFILE,"# shift end time [int shift_end]\n");
   fprintf(outFILE,"%d\n",modelData->shift_end);

   /* int*                  t_service;              /\**< nC-array with (average) customer service times (in seconds) *\/ */
   fprintf(outFILE,"# service times at customers (depot has service time 0.0) [int* t_service]\n");
   assert( modelData->t_service != NULL );
   for( i = 0; i < modelData->nC; i++ )
   {
      fprintf(outFILE,"%d ",modelData->t_service[i]);
   }
   fprintf(outFILE,"\n");

   /* int*                  t_service_maxDev;       /\**< nC-array with estimated maximal deviation of customer service times beyond average (in seconds) *\/ */
   fprintf(outFILE,"# max. deviation of service times at customers (depot again 0.0) [int* t_service_maxDev]\n");
   assert( modelData->t_service_maxDev != NULL );
   for( i = 0; i < modelData->nC; i++ )
   {
      fprintf(outFILE,"%d ",modelData->t_service_maxDev[i]);
   }
   fprintf(outFILE,"\n");

   /* int**                 t_travelAM;             /\**< (nC x nC)-array with morning's travel times between customers in data set (in seconds) *\/ */
   fprintf(outFILE,"# morning/AM travel time matrix between customers (incl. depot as last customer) [int** t_travelAM]\n");
   assert( modelData->t_travelAM != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      assert( modelData->t_travelAM[i] != NULL );
      for( j = 0; j < modelData->nC; j++ ) /* ...to */
      {
         fprintf(outFILE,"%d ",modelData->t_travelAM[i][j]);
      }
      fprintf(outFILE,"\n");
   }

   /* int**                 t_travelNoon;           /\**< (nC x nC)-array with noon's travel times between customers in data set (in seconds) *\/ */
   fprintf(outFILE,"# noon travel time matrix between customers (incl. depot as last customer) [int** t_travelNoon]\n");
   assert( modelData->t_travelNoon != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      assert( modelData->t_travelNoon[i] != NULL );
      for( j = 0; j < modelData->nC; j++ ) /* ...to */
      {
         fprintf(outFILE,"%d ",modelData->t_travelNoon[i][j]);
      }
      fprintf(outFILE,"\n");
   }

   /* int**                 t_travelPM;             /\**< (nC x nC)-array with afternoon's travel times between customers in data set (in seconds) *\/ */
   fprintf(outFILE,"# afternoon/PM travel time matrix between customers (incl. depot as last customer) [int** t_travelPM]\n");
   assert( modelData->t_travelPM != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      assert( modelData->t_travelPM[i] != NULL );
      for( j = 0; j < modelData->nC; j++ ) /* ...to */
      {
         fprintf(outFILE,"%d ",modelData->t_travelPM[i][j]);
      }
      fprintf(outFILE,"\n");
   }

   /* int**                 t_travel;               /\**< (nC x nC)-array with average travel times between customers in data set (in seconds) *\/ */
   fprintf(outFILE,"# average travel time matrix between customers (incl. depot as last customer) [int** t_travel]\n");
   assert( modelData->t_travel != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      assert( modelData->t_travel[i] != NULL );
      for( j = 0; j < modelData->nC; j++ ) /* ...to */
      {
         fprintf(outFILE,"%d ",modelData->t_travel[i][j]);
      }
      fprintf(outFILE,"\n");
   }

   /* int**                 t_travel_maxDev;        /\**< (nC x nC)-array with estimated maximal deviation of travel times between customers beyond average (in seconds) *\/ */
   fprintf(outFILE,"# max. deviation of (average) travel time matrix between customers (incl. depot as last customer) [int** t_travel_maxDev]\n");
   assert( modelData->t_travel_maxDev != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      assert( modelData->t_travel_maxDev[i] != NULL );
      for( j = 0; j < modelData->nC; j++ ) /* ...to */
      {
         fprintf(outFILE,"%d ",modelData->t_travel_maxDev[i][j]);
      }
      fprintf(outFILE,"\n");
   }

   /* SCIP_Bool             workOnSaturdays;        /\**< indicator whether saturdays are to be counted as workdays or not (TRUE/FALSE (1/0 unsigned int.)); for information/output purposes *\/ */
   fprintf(outFILE,"# indicator whether saturdays count as workdays (stored as 1 or 0 for \"TRUE\" or \"FALSE\", resp.) [SCIP_Bool workOnSaturdays]\n");
   fprintf(outFILE,"%d\n",(int)modelData->workOnSaturdays);

   /* date*                 startDate;              /\**< struct containing the start date of the planning period (day, month, year, weekday); for information/output purposes *\/ */
   fprintf(outFILE,"# start date of planning period (stored as 4 integers: day, month, year, weekday) [date* startDate]\n");
   assert( modelData->startDate != NULL );
   fprintf(outFILE,"%d %d %d %d\n",(modelData->startDate)->day,(modelData->startDate)->month,(modelData->startDate)->year,(modelData->startDate)->weekday);

   /* date*                 endDate;                /\**< struct containing the end date of the planning period (day, month, year, weekday); for information/output purposes *\/ */
   fprintf(outFILE,"# end date of planning period (stored as 4 integers: day, month, year, weekday) [date* endDate]\n");
   assert( modelData->endDate != NULL );
   fprintf(outFILE,"%d %d %d %d\n",(modelData->endDate)->day,(modelData->endDate)->month,(modelData->endDate)->year,(modelData->endDate)->weekday);

   /* int                   maxDelayEvents;         /\**< maximal number of delay events (deviations in service and travel times combined) to be "robustified" against *\/ */
   fprintf(outFILE,"# maximal number of delay events considered (in robust optimization) [int maxDelayEvents]\n");
   fprintf(outFILE,"%d\n",modelData->maxDelayEvents);

   /* int                   nDays;                  /\**< number of work days in planning period (from startDate to endDate, including both) *\/ */
   fprintf(outFILE,"# number of workdays in planning period [int nDays]\n");
   fprintf(outFILE,"%d\n",modelData->nDays);

   /* dayProperties*        days;                   /\**< struct containing the properties of each day *\/ */
   fprintf(outFILE,"# day properties (stored as 6 integers: index, shift_start, shift_end, break_start, break_end, break_length) [dayProperty* days]\n");
   for( i = 0; i < modelData->nDays; i++)
   {
       fprintf(outFILE, "%d ",modelData->days[i]->index);
       fprintf(outFILE, "%d ",modelData->days[i]->shift_start);
       fprintf(outFILE, "%d ",modelData->days[i]->shift_end);
       fprintf(outFILE, "%d ",modelData->days[i]->break_start);
       fprintf(outFILE, "%d ",modelData->days[i]->break_end);
       fprintf(outFILE, "%d\n",modelData->days[i]->break_length);
   }
    fprintf(outFILE,"# day sizes (day, number of available customers on that day) [int* day_sizes]\n");
   for( i = 0; i < modelData->nDays; i++)
   {
       fprintf(outFILE, "%d ",i);
       fprintf(outFILE, "%d\n",modelData->day_sizes[i]);
   }

   /* neighbor***           neighbors;              /\**< (nC x nDays) array of linked lists with neighbor nodes of customer nodes in total network (can be NULL <-> complete digraph) *\/ */
   fprintf(outFILE,"# list of neighbor nodes of customers for each day in planning period [neighbor*** neighbors]\n");
   /* rows are printed for each customer only for the days on which that customer has a positive number of neighbors: customer index, (model) day index, neighbor list */
   if( modelData->neighbors == NULL )
   {
      fprintf(outFILE,"-1\n"); /* to indicate neighbor list is not present <-> complete digraph for every day */
   }
   else
   {
      assert( modelData->neighbors != NULL );
      for( i = 0; i < modelData->nC; i++ )
      {
         assert( modelData->neighbors[i] != NULL );
         for( j = 0; j < modelData->nDays; j++ )
         {
            if( modelData->neighbors[i][j] != NULL ) /* no output for "empty days" */
            {
               tmpnb = modelData->neighbors[i][j];
               fprintf(outFILE,"%d %d ",i,j);
               while( tmpnb != NULL )
               {
                  fprintf(outFILE,"%d ",tmpnb->id);
                  tmpnb = tmpnb->next;
               }
               fprintf(outFILE,"\n");
            }
         }
      }
   }

   fclose(outFILE);

   return SCIP_OKAY;
}



/* read model data from file */
SCIP_RETCODE readModelData(
   SCIP* scip,                            /**< SCIP pointer */
   const char* inputFile,                 /**< input file name */
   model_data* modelData                  /**< pointer to model data struct (alloc'd outside of this function!) */
   )
{
   FILE* inFILE = NULL;
   int i,j,k,c;
   modelWindow* tmptw = NULL;
   modelWindow* twHead = NULL;
   int MAXSTRLEN = 100*SCIP_MAXSTRLEN; /* need more than SCIP_MAXSTRLEN = 1024 at least for customer ID line */
   char strbuffer[MAXSTRLEN];
   char* token = NULL;
   neighbor* tmpnb = NULL;
   neighbor* nbHead = NULL;

   assert( modelData != NULL );
   assert( inputFile != NULL );
   inFILE = fopen(inputFile,"r");
   assert( inFILE != NULL );

   /* format being read is assumed to be EXACTLY as produced by writeModelData -> only rudimentary checks are performed to ensure compatability */

   /* nC (number of customers + 1) */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int nC */
   assert( ! strncmp(strbuffer,"# number of customers",21) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&(modelData->nC));
   assert( k == 1 );

   /* customerIDs */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int* customerIDs */
   assert( ! strncmp(strbuffer,"# customer IDs",14) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->customerIDs),modelData->nC) );
   assert( modelData->customerIDs != NULL );
   k = 0;
   i = 0;
   token = strtok(strbuffer," ");
   while( token != NULL )
   {
      k += sscanf(token,"%d",&(modelData->customerIDs[i++]));
      token = strtok(NULL," ");
      if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
      {
         break;
      }
   }
   assert( k == modelData->nC );
   assert( i == k );

   /* time windows */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field modelWindows** timeWindows */
   assert( ! strncmp(strbuffer,"# time window lists",19) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->timeWindows),modelData->nC) );
   assert( modelData->timeWindows != NULL );
   for( i = 0; i < modelData->nC-1; i++ ) /* depot gets NULL */
   {
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      token = strtok(strbuffer," ");
      k = sscanf(token,"%d",&j);
      assert( k == 1 );
      k = 0;
      c = 0;
      token = strtok(NULL," ");
      /* assert( j > 0 ); */ /* use only if assumption is that every customer has at least one time window available for
         being serviced; see corresp. comment in writeModelData */
      if( j > 0 )
      {
         twHead = NULL;
         while( token != NULL )
         {
            assert( tmptw == NULL );
            SCIP_CALL( SCIPallocBlockMemory(scip,&tmptw) );
            assert( tmptw != NULL );
            k += sscanf(token,"%d",&(tmptw->day));
            token = strtok(NULL," ");
            k += sscanf(token,"%d",&(tmptw->start_t));
            token = strtok(NULL," ");
            k += sscanf(token,"%d",&(tmptw->end_t));
            token = strtok(NULL," ");
            k += sscanf(token,"%lf",&(tmptw->weigth));
            token = strtok(NULL," ");
            assert( c == k-4 );
            c = k;
            tmptw->next = NULL;
            if( c == 4 ) /* first list element was just constructed */
            {
               (modelData->timeWindows[i]) = tmptw;
               twHead = tmptw;
            }
            else /* append to list */
            {
               twHead->next = tmptw;
               twHead = twHead->next;
            }
            assert( twHead != NULL );
            tmptw = tmptw->next;
            if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
            {
               break;
            }
         }
         assert( c == 4*j );
      }
      else /* empty time window list for this customer... should probably not happen, see earlier comments */
      {
         assert( ! strcmp(token,"\n") ); /* only a trailing space is left in this line */
         modelData->timeWindows[i] = NULL;
      }
   }
   /* depot */
   SCIP_CALL( SCIPallocBlockMemory(scip,&tmptw) );
   tmptw->day = 0;
   tmptw->start_t = DEFAULT_SHIFT_START;
   tmptw->end_t = DEFAULT_SHIFT_END;
   tmptw->next = NULL;
   tmptw->weigth = 1.0;
   modelData->timeWindows[modelData->nC - 1] = tmptw;
   twHead = tmptw;


   /* objective coefficients */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field SCIP_Real* obj */
   assert( ! strncmp(strbuffer,"# objective",11) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->obj),modelData->nC) );
   assert( modelData->obj != NULL );
   k = 0;
   i = 0;
   token = strtok(strbuffer," ");
   while( token != NULL )
   {
      k += sscanf(token,"%lf",&(modelData->obj[i++]));
      token = strtok(NULL," ");
      if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
      {
         break;
      }
   }
   assert( k == modelData->nC );
   assert( i == k );

   /* shift start time */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int shift_start */
   assert( ! strncmp(strbuffer,"# shift start",13) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&(modelData->shift_start));
   assert( k == 1 );

   /* shift end time */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int shift_end */
   assert( ! strncmp(strbuffer,"# shift end",11) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&(modelData->shift_end));
   assert( k == 1 );

   /* (average) service times */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int* t_service */
   assert( ! strncmp(strbuffer,"# service times",15) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_service),modelData->nC) );
   assert( modelData->t_service != NULL );
   k = 0;
   i = 0;
   token = strtok(strbuffer," ");
   while( token != NULL )
   {
      k += sscanf(token,"%d",&(modelData->t_service[i++]));
      token = strtok(NULL," ");
      if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
      {
         break;
      }
   }
   assert( k == modelData->nC );
   assert( i == k );

   /* estim. max. deviation of service times */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int* t_service_maxDev */
   assert( ! strncmp(strbuffer,"# max. deviation of service times",33) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_service_maxDev),modelData->nC) );
   assert( modelData->t_service_maxDev != NULL );
   k = 0;
   i = 0;
   token = strtok(strbuffer," ");
   while( token != NULL )
   {
      k += sscanf(token,"%d",&(modelData->t_service_maxDev[i++]));
      token = strtok(NULL," ");
      if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
      {
         break;
      }
   }
   assert( k == modelData->nC );
   assert( i == k );

   /* morning travel time matrix */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int** t_travelAM */
   assert( ! strncmp(strbuffer,"# morning/AM travel time",24) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelAM),modelData->nC) );
   assert( modelData->t_travelAM != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelAM[i]),modelData->nC) );
      assert( modelData->t_travelAM[i] != NULL );
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      k = 0;
      j = 0;
      token = strtok(strbuffer," ");
      while( token != NULL ) /* ... to */
      {
         k += sscanf(token,"%d",&(modelData->t_travelAM[i][j++]));
         token = strtok(NULL," ");
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
         {
            break;
         }
      }
      assert( k == modelData->nC );
      assert( j == k );
   }

   /* noon travel time matrix */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int** t_travelNoon */
   assert( ! strncmp(strbuffer,"# noon travel time",18) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelNoon),modelData->nC) );
   assert( modelData->t_travelNoon != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelNoon[i]),modelData->nC) );
      assert( modelData->t_travelNoon[i] != NULL );
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      k = 0;
      j = 0;
      token = strtok(strbuffer," ");
      while( token != NULL ) /* ... to */
      {
         k += sscanf(token,"%d",&(modelData->t_travelNoon[i][j++]));
         token = strtok(NULL," ");
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
         {
            break;
         }
      }
      assert( k == modelData->nC );
      assert( j == k );
   }

   /* afternoon travel time matrix */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int** t_travelPM */
   assert( ! strncmp(strbuffer,"# afternoon/PM travel time",26) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelPM),modelData->nC) );
   assert( modelData->t_travelPM != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travelPM[i]),modelData->nC) );
      assert( modelData->t_travelPM[i] != NULL );
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      k = 0;
      j = 0;
      token = strtok(strbuffer," ");
      while( token != NULL ) /* ... to */
      {
         k += sscanf(token,"%d",&(modelData->t_travelPM[i][j++]));
         token = strtok(NULL," ");
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
         {
            break;
         }
      }
      assert( k == modelData->nC );
      assert( j == k );
   }

   /* average travel time matrix */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int** t_travel */
   assert( ! strncmp(strbuffer,"# average travel time",21) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travel),modelData->nC) );
   assert( modelData->t_travel != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travel[i]),modelData->nC) );
      assert( modelData->t_travel[i] != NULL );
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      k = 0;
      j = 0;
      token = strtok(strbuffer," ");
      while( token != NULL ) /* ... to */
      {
         k += sscanf(token,"%d",&(modelData->t_travel[i][j++]));
         token = strtok(NULL," ");
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
         {
            break;
         }
      }
      assert( k == modelData->nC );
      assert( j == k );
   }

   /* max. deviations of travel times */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int* t_travel_maxDev */
   assert( ! strncmp(strbuffer,"# max. deviation of (average) travel time",40) );
   SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travel_maxDev),modelData->nC) );
   assert( modelData->t_travel_maxDev != NULL );
   for( i = 0; i < modelData->nC; i++ ) /* from... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->t_travel_maxDev[i]),modelData->nC) );
      assert( modelData->t_travel_maxDev[i] != NULL );
      (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
      k = 0;
      j = 0;
      token = strtok(strbuffer," ");
      while( token != NULL ) /* ... to */
      {
         k += sscanf(token,"%d",&(modelData->t_travel_maxDev[i][j++]));
         token = strtok(NULL," ");
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak */
         {
            break;
         }
      }
      assert( k == modelData->nC );
      assert( j == k );
   }

   /* workOnSaturdays indicator */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field SCIP_Bool workOnSaturdays */
   assert( ! strncmp(strbuffer,"# indicator whether saturdays",29) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&j);
   assert( k == 1 );
   assert( (j == 0) || (j == 1) );
   modelData->workOnSaturdays = (SCIP_Bool)j;

   /* start date */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field date* startDate */
   assert( ! strncmp(strbuffer,"# start date",12) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemory(scip,&(modelData->startDate)) );
   assert( modelData->startDate != NULL );
   k = sscanf(strbuffer,"%d %d %d %d",&((modelData->startDate)->day),&((modelData->startDate)->month),&((modelData->startDate)->year),&((modelData->startDate)->weekday));
   assert( k == 4 );

   /* end date */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field date* endDate */
   assert( ! strncmp(strbuffer,"# end date",10) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   SCIP_CALL( SCIPallocBlockMemory(scip,&(modelData->endDate)) );
   assert( modelData->endDate != NULL );
   k = sscanf(strbuffer,"%d %d %d %d",&((modelData->endDate)->day),&((modelData->endDate)->month),&((modelData->endDate)->year),&((modelData->endDate)->weekday));
   assert( k == 4 );

   /* max. number of delay events */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int maxDelayEvents */
   assert( ! strncmp(strbuffer,"# maximal number of delay events",32) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&(modelData->maxDelayEvents));
   assert( k == 1 );

   /* number of workdays in planning period */
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field int nDays */
   assert( ! strncmp(strbuffer,"# number of workdays",20) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   k = sscanf(strbuffer,"%d",&(modelData->nDays));
   assert( k == 1 );

   /* time windows for depot */
    for( i = 1; i < modelData->nDays; i++ )
    {
        SCIP_CALL( SCIPallocBlockMemory(scip,&tmptw) );
        tmptw->day = i;
        tmptw->start_t = DEFAULT_SHIFT_START;
        tmptw->end_t = DEFAULT_SHIFT_END;
        tmptw->weigth = 1.0;
        tmptw->next = NULL;
        twHead->next = tmptw;
        twHead = twHead->next;
    }

    /* day properties */
    SCIP_CALL( SCIPallocBlockMemoryArray( scip, &(modelData->days), modelData->nDays) );
    (void)! fgets(strbuffer, MAXSTRLEN, inFILE); /* commentary line -> field dayProperty* days */
    if(! strncmp(strbuffer,"# day properties",16) )
    {
        for ( i = 0; i < modelData->nDays; i++)
        {
            SCIP_CALL( SCIPallocBlockMemory(scip, &(modelData->days[i])) );
            k = 0;
            (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
            token = strtok(strbuffer, " ");
            k += sscanf(token, "%d", &j);
            assert(j == i);
            modelData->days[i]->index = i;
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->days[i]->shift_start));
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->days[i]->shift_end));
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->days[i]->break_start));
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->days[i]->break_end));
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->days[i]->break_length));
            assert(k == 6);
        }
        (void)! fgets(strbuffer, MAXSTRLEN, inFILE);/* commentary line -> field neighbor*** neighbors */
    }
    /* day sizes */
    SCIP_CALL( SCIPallocBlockMemoryArray( scip, &(modelData->day_sizes), modelData->nDays) );
    if( ! strncmp(strbuffer,"# day sizes",11) )
    {
        for ( i = 0; i < modelData->nDays; i++)
        {
            k = 0;
            (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
            token = strtok(strbuffer, " ");
            k += sscanf(token, "%d", &j);
            assert(j == i);
            token = strtok(NULL, " ");
            k += sscanf(token, "%d", &(modelData->day_sizes[i]));
            assert(k == 2);
        }
        (void)! fgets(strbuffer, MAXSTRLEN, inFILE);/* commentary line -> field neighbor*** neighbors */
    }else
    {
        for ( i = 0; i < modelData->nDays; i++)
        {
            modelData->day_sizes[i] = modelData->nC - 1;
        }
    }
   /* list of neighbor nodes of customers for each day in planning period */
   assert( ! strncmp(strbuffer,"# list of neighbor nodes",24) );
   (void)! fgets(strbuffer, MAXSTRLEN, inFILE);
   /* -1 if no neighbor lists are given; else rows contain: customer index, (model) day index, neighbor list */
   if( ! strncmp(strbuffer,"-1",2) )
   {
      modelData->neighbors = NULL;
   }
   else /* at least one neighbor list is given -> allocate and initialize list array, read the first one... */
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->neighbors),modelData->nC) );
      assert( modelData->neighbors != NULL );
      for( i = 0; i < modelData->nC; i++ )
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip,&(modelData->neighbors[i]),modelData->nDays) );
         assert( modelData->neighbors[i] != NULL );
         for( j = 0; j < modelData->nDays; j++ )
         {
            modelData->neighbors[i][j] = NULL; /* init. empty lists */
         }
      }
      k = 0;
      token = strtok(strbuffer," ");
      k = sscanf(token,"%d",&i); /* customer index (start node) */
      assert( k == 1 );
      token = strtok(NULL," ");
      k = sscanf(token,"%d",&j); /* day index */
      assert( k == 1 );
      token = strtok(NULL," ");
      assert( strcmp(token,"\n") != 0 );
      assert( (i >= 0) && (i < modelData->nC) && (j >= 0) && (j < modelData->nDays) );
      SCIP_CALL( SCIPallocBlockMemory(scip,&tmpnb) );
      assert( tmpnb != NULL );
      nbHead = modelData->neighbors[i][j];
      assert( nbHead == NULL );
      while( token != NULL )
      {
         k = sscanf(token,"%d",&(tmpnb->id));
         assert( k == 1 );
         token = strtok(NULL," ");
         tmpnb->next = NULL;
         if( nbHead == NULL ) /* first list element */
         {
            nbHead = tmpnb;
            modelData->neighbors[i][j] = tmpnb;
            tmpnb = tmpnb->next;
            assert( nbHead != NULL );
            assert( tmpnb == NULL );
         }
         else /* append to list */
         {
            nbHead->next = tmpnb;
            nbHead = nbHead->next;
            assert( nbHead != NULL );
            tmpnb = tmpnb->next;
            assert( tmpnb == NULL );
         }
         if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak (reached end of list) */
         {
            break;
         }
         /* else: more elements are coming up */
         SCIP_CALL( SCIPallocBlockMemory(scip,&tmpnb) );
         assert( tmpnb != NULL );
      }
      /* now read rest of the neighbor lists (until end of file is reached) */
      while( fgets(strbuffer, MAXSTRLEN, inFILE) ) /* stops at EOF */
      {
         k = 0;
         token = strtok(strbuffer," ");
         k = sscanf(token,"%d",&i); /* customer index (start node) */
         assert( k == 1 );
         token = strtok(NULL," ");
         k = sscanf(token,"%d",&j); /* day index */
         assert( k == 1 );
         token = strtok(NULL," ");
         assert( strcmp(token,"\n") != 0 );
         assert( (i >= 0) && (i < modelData->nC) && (j >= 0) && (j < modelData->nDays) );
         SCIP_CALL( SCIPallocBlockMemory(scip,&tmpnb) );
         assert( tmpnb != NULL );
         nbHead = modelData->neighbors[i][j];
         assert( nbHead == NULL );
         while( token != NULL )
         {
            k = sscanf(token,"%d",&(tmpnb->id));
            assert( k == 1 );
            token = strtok(NULL," ");
            tmpnb->next = NULL;
            if( nbHead == NULL ) /* first list element */
            {
               nbHead = tmpnb;
               modelData->neighbors[i][j] = tmpnb;
               tmpnb = tmpnb->next;
               assert( nbHead != NULL );
               assert( tmpnb == NULL );
            }
            else /* append to list */
            {
               nbHead->next = tmpnb;
               nbHead = nbHead->next;
               assert( nbHead != NULL );
               tmpnb = tmpnb->next;
               assert( tmpnb == NULL );
            }
            if( ! strcmp(token,"\n") ) /* catch trailing space -> linebreak (reached end of list) */
            {
               break;
            }
            /* else: more elements are coming up */
            SCIP_CALL( SCIPallocBlockMemory(scip,&tmpnb) );
            assert( tmpnb != NULL );
         }
      }
   }



   fclose(inFILE);

   return SCIP_OKAY;
}
