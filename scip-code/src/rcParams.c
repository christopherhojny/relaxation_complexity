/**@file   rcParams.cpp
 * @brief  set up parameters for computing RC
 * @author Christopher Hojny
 */

#include "rcParams.h"


/* default settings for parameters - for documentation see below */

/* information regarding problem */
#define DEFAULT_PRINTDATAINFO                FALSE
#define DEFAULT_HANDLESYMMETRY                TRUE
#define DEFAULT_HIDINGSETCUTS                 TRUE
#define DEFAULT_METHOD                           2
#define DEFAULT_EPSILON                      0.001

/** Set basic SCIP parameters that are relevant for computing RC */
extern
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   return SCIP_OKAY;
}


/** Introduce parameters that are relevant for computing RC */
extern
SCIP_RETCODE addRCParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   /* information regarding problem */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/printdatainfo", "Shall information on the data points be printed?",
         NULL, TRUE, DEFAULT_PRINTDATAINFO, NULL, NULL) );

   /* type of method to solve problem */
   SCIP_CALL( SCIPaddIntParam(scip, "rc/method", "Which method shall be used to find RC (0: compact MIP, 1: column generation, 2: hybrid)",
         NULL, TRUE, DEFAULT_METHOD, 0, 2, NULL, NULL) );

   /* whether we shall use symmetry handling techniques */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/handlesymmetry", "Shall symmetry be handled in compact model?",
         NULL, TRUE, DEFAULT_HANDLESYMMETRY, NULL, NULL) );

   /* whether we shall use cuts based on hiding sets */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/hidingsetcuts", "Shall hiding set based cuts be added to the model?",
         NULL, TRUE, DEFAULT_HIDINGSETCUTS, NULL, NULL) );

   /* the epsilon used to rate an inequality as violated */
   SCIP_CALL( SCIPaddRealParam(scip, "rc/epsilon", "The epsilon used to rate an inequality as violated",
         NULL, TRUE, DEFAULT_EPSILON, 0.00001, 1.0, NULL, NULL) );

   return SCIP_OKAY;
}
