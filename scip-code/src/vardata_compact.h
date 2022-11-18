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

/**@file   vardata_compact.h
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_COMPACT__
#define __SCIP_VARDATA_COMPACT__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates is violated variable */
SCIP_RETCODE SCIPcreateVarIsviolatedCompact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

/** create variable data */
SCIP_RETCODE SCIPvardataCreateCompact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int                   idx                 /**< index of corresponding inequality */
   );

/** get index of corresponding inequality */
int SCIPvardataGetInequalityidx(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

#ifdef __cplusplus
}
#endif

#endif
