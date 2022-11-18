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

/**@file   rcParams.h
 * @brief  set up parameters for computing RC
 * @author Christopher Hojny
 */

#ifndef __RCPARAMS_H__
#define __RCPARAMS_H__

// SCIP include
#include <scip/scip.h>

#ifdef __cplusplus
extern "C" {
#endif

/** Set basic SCIP parameters that are relevant for computing RC */
extern
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** Introduce parameters that are relevant for computing RC */
extern
SCIP_RETCODE addRCParameters(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
