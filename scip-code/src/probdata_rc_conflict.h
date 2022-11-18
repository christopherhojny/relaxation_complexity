/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program computeRC                             */
/*                                                                           */
/*    an implementation of a branch-and-cut and branch-and-price             */
/*    algorithm to compute the epsilon relaxation complexity of              */
/*    a full-dimensional lattice-convex set X and a finite set               */
/*    of points Y.                                                           */
/*                                                                           */
/*    Copyright (C) 2022-     Gennadiy Averkov, Christopher Hojny,           */
/*                            Matthias Schymura                              */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/*    Copyright (C) 2002-2022 Zuse Institute Berlin                          */
/*                                                                           */
/*       mailto: scip@zib.de                                                 */
/*       Licensed under the Apache License, Version 2.0                      */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

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
   int                   absmax,             /**< maximum absolute value of a coordinate in X */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities (or NULL if not needed) */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality
                                                (or NULL if not needed) */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality
                                                (or NULL if not needed) */
   int                   ninequalities,      /**< number of inequalities encoded in previous data structures
                                                (or -1 if not needed) */
   int                   maxninequalities    /**< maximum number of inequalities that can be stored
                                                (if allocated) */
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
