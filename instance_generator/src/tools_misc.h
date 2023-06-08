/**@file   tools_misc.h
 * @brief  miscellaneous auxiliary functions
 * @author Andreas M. Tillmann, TU Braunschweig
 */

#ifndef __TOOLS_MISC_H__
#define __TOOLS_MISC_H__

#include <scip/scip.h>



/** read comand line arguments */
extern
SCIP_RETCODE readArguments(int argc, char** argv, char** baseDataPath, char** dataSetName, char** customerFile, char** matrequestFile, char** matresponseAMFile, char** matresponseNoonFile, char** matresponsePMFile, char** metaFile, char** outputFile);



#endif
