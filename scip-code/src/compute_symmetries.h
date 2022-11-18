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

/**@file   compute_symmetries.h
 * @brief  functions to compute symmetries
 * @author Christopher Hojny
 */

#ifndef __COMPUTE_SYMMETRIES_H__
#define __COMPUTE_SYMMETRIES_H__

#include "datapoints.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


extern
/** computes common symmetries of X and Y */
SCIP_RETCODE computeSymmetries(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int***                perms,              /**< pointer to store permutations */
   int*                  nperms,             /**< pointer to store number of permutations stored in perms */
   int*                  nmaxperms           /**< pointer to store maximum number of permutations fiiting in perms */
   );

#ifdef __cplusplus
}
#endif

#endif
