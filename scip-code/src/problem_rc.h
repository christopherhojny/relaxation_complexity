/**@file   problem_clustering.h
 * @brief  Problem data for computing RC
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __PROBLEM_CLUSTERING_H__
#define __PROBLEM_CLUSTERING_H__

#include "datapoints.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates initial model for computing RC */
extern
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmax,             /**< maximum absolute value of a coordinate in X */
   SCIP_Bool             secondphase         /**< whether we are in the second phase of the hybrid approach */
   );

/** free problem data */
extern
SCIP_RETCODE SCIPfreeModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y                   /**< pointer to data points of Y */
   );

/** solves the problem */
extern
SCIP_RETCODE SCIPsolveRCproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   ub,                 /**< upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmax,             /**< maximum absolute value of a coordinate in X */
   char*                 filename,           /**< name of file for storing solution (or NULL) */
   SCIP_Real             timelimit           /**< time limit */
   );

#ifdef __cplusplus
}
#endif

#endif

