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

/**@file   typedefs.h
 * @brief  macros used throughout the entire project
 * @author Christopher Hojny
 */

#ifndef __TYPEDEFS_H_
#define __TYPEDEFS_H_

#ifdef __cplusplus
extern "C" {
#endif

#define METHOD_COMPACT_MIP      0                   /**< compact MIP model */
#define METHOD_CG               1                   /**< column generation */
#define METHOD_HYBRID_COMPACT   2                   /**< column generation for generating bound and compact for rest */
#define METHOD_CONFLICT         3                   /**< conflict based model */
#define METHOD_HYBRID_CONFLICT  4                   /**< column generation for generating bound and conflict for rest */

#define GREEDY_SORT_STANDARD    0                   /**< standard ordering of infeasible points for greedy heuristic */
#define GREEDY_SORT_DISTANCE    1                   /**< ordering of infeasible points based on lattice distance to feasible points */

#ifdef __cplusplus
}
#endif

#endif
