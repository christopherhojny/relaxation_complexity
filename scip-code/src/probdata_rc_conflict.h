/**@file   probdata_rc_conflict.h
 * @brief  Problem data for computing RC using a conflict based model
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROBDATA_RC_CONFLICT_H__
#define __PROBDATA_RC_CONFLICT_H__

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmax              /**< maximum absolute value of a coordinate in X */
   );

/**< returns feasible points */
Datapoints* SCIPprobdataConflictGetX(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/**< returns infeasible points */
Datapoints* SCIPprobdataConflictGetY(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/**< returns upper bound on number of classes */
int SCIPprobdataConflictGetUb(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/**< returns violation variables */
SCIP_VAR*** SCIPprobdataConflictGetViolatedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

/**< returns isused variables */
SCIP_VAR** SCIPprobdataConflictGetIsusedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   );

#ifdef __cplusplus
}
#endif

#endif
