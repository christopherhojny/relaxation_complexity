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
 * @brief  file implementing functions for convex hull computations
 * @author Christopher Hojny
 */

#ifndef __CONVEX_HULL_H__
#define __CONVEX_HULL_H__

#include "datapoints.h"
#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/* constructs a CDD matrix containing points from a set of data points */
extern
dd_MatrixPtr constructGeneratorMatrixPoints(
   Datapoints*           datapoints,         /**< data points */
   int*                  selectedpoints,     /**< array of selected data points */
   int                   nselectedpoints     /**< number of selected data points */
   );

/** computes facet description of set of points */
extern
dd_MatrixPtr computeConvexHullFacets(
   SCIP*                 scip,               /**< SCIP instance */
   dd_MatrixPtr          generators,         /**< points generating convex hull */
   SCIP_Bool*            success             /**< pointer to store whether we were successful */
   );

#ifdef __cplusplus
}
#endif

#endif
