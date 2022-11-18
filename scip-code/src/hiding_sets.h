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
 * @brief  functions for hiding set computations
 * @author Christopher Hojny
 */

#ifndef __HIDING_SETS_H__
#define __HIDING_SETS_H__

#include "datapoints.h"
#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


extern
/** computes index sets of hiding set cuts, i.e., pairs of integer points forming a hiding set */
SCIP_RETCODE computeHidingSetCuts(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int**                 hidingsetidx,       /**< address to pointer for storing hiding set indices, i.e.,
                                              *   (tuples) of indices of infeasible points */
   int*                  nhidingsetidx,      /**< number of hiding set indices encoded in hidingsetidx */
   int*                  maxnhidingsetidx    /**< maximum number of entries in hidingsetidx */
   );

#ifdef __cplusplus
}
#endif

#endif
