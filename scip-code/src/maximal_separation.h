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

/**@file   convex_hull.h
 * @brief  functions for computing maximally separable sets
 * @author Christopher Hojny
 */

#ifndef __MAXIMAL_SEPARATION_H__
#define __MAXIMAL_SEPARATION_H__

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** finds an inequality separating as many infeasible points as possible */
SCIP_EXPORT
SCIP_RETCODE SCIPfindMaximalSeparatingInequality(
   SCIP*                 scip,               /**< main SCIP pointer */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   nX,                 /**< number of data points encoded in X */
   int                   nY,                 /**< number of data points encoded in Y */
   SCIP_Real*            objvals,            /**< objective values of points to be separated (or NULL for all ones) */
   int                   dimension,          /**< dimension of data points */
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   int*                  fixedY,             /**< array of points in Y whose separation status is fixed
                                                  (or NULL if not needed)*/
   SCIP_Real*            statusfixedY,       /**< separation status of points in fixedY (or NULL if not needed):
                                                0 - not separated, 1 - separated */
   int                   nfixedY,            /**< number of points in fixedY */
   int*                  separatedY,         /**< array of points separated in solution (needs to be allocated) */
   int*                  nseparatedY,        /**< pointer to store number of separated points */
   SCIP_Real*            sepainequality,     /**< array to store separated inequality (needs to be allocated) */
   SCIP_CONSHDLR*        conshdlr,           /**< conshdlr for branching decisions (or NULL if not needed) */
   SCIP_Longint          nodelimit,          /**< node limit for solving the separation problem (or -1 if no limit) */
   SCIP_Real             timelimit,          /**< time limit for finding maximumally violated inequality (of -1 if no limit) */
   SCIP_Bool*            success             /**< pointer to store whether a separating inequality could be found */
   );

/** finds a greedy solution for RC */
SCIP_EXPORT
SCIP_RETCODE heurGreedy(
   SCIP*                 scip,               /**< SCIP pointer */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   nX,                 /**< number of data points encoded in X */
   int                   nY,                 /**< number of data points encoded in Y */
   int                   dimension,          /**< dimension of data points */
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   int                   ub,                 /**< upper bound on RC */
   int                   greedysort,         /**< sorting of infeasible points used by heuristic */
   int                   method,             /**< method used to compute RC */
   SCIP_Real**           inequalities,       /**< (ub x dimension+1) array to store inequalities */
   int**                 separatedpoints,    /**< (ub x nY) array to store list of separated points */
   int*                  nseparatedpoints,   /**< ub-dimensional array to store for each inequality
                                                the number of separated points */
   int*                  ninequalities,      /**< pointer to store number of inequalities in solution */
   dd_MatrixPtr          facetsconvexhull,   /**< facets of conv(X) (or NULL if not needed) */
   SCIP_Real             timelimit           /**< time limit for greedy computation */
   );


#ifdef __cplusplus
}
#endif

#endif
