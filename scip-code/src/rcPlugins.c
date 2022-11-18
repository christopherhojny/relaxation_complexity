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

/**@file   rcPlugins.c
 * @brief  load SCIP plugins for computing RC
 * @author Christopher Hojny
 */

#include "cons_conflict.h"
#include "rcPlugins.h"
#include "prop_convexity.h"
#include "prop_intersection.h"

#include "scip/scipdefplugins.h"


/** Include basic plugins needed for computing RC */
SCIP_RETCODE includeRCPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   SCIP_CALL( SCIPincludePropConvexity(scip) );
   SCIP_CALL( SCIPincludePropIntersection(scip) );
   SCIP_CALL( SCIPincludeConshdlrConflict(scip) );

   return SCIP_OKAY;
}
