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

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"

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
   SCIP_Bool             secondphase,        /**< whether we are in the second phase of the hybrid approach */
   SCIP_Real*            heurtime,           /**< pointer to store time spent in heuristics */
   SCIP_Bool*            success             /**< pointer to store whether model could be built successfully */
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

/** returns solution of the problem based on conv(X) based on CDD representation */
extern
SCIP_RETCODE SCIPgetSolutionCompactModelCDD(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   dd_MatrixPtr          facetsconvexhull    /**< facet description of conv(X) */
   );

/** fills existing solution of the problem based on explicit solution */
extern
SCIP_RETCODE SCIPgetSolutionCompactModelExplicit(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality */
   int                   ninequalities       /**< number of inequalities encoded in previous data structures */
   );

/** fills existing solution of the problem based on conv(X) based on CDD representation */
extern
SCIP_RETCODE SCIPgetSolutionConflictModelCDD(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   dd_MatrixPtr          facetsconvexhull    /**< facet description of conv(X) */
   );

/** fills existing solution of the problem based on explicit solution */
extern
SCIP_RETCODE SCIPgetSolutionConflictModelExplicit(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality */
   int                   ninequalities       /**< number of inequalities encoded in previous data structures */
   );

#ifdef __cplusplus
}
#endif

#endif

