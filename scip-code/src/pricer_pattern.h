/*
  This file is based on the file pricer_binpack.h distributed by via SCIP:

  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  *                                                                           *
  *                  This file is part of the program and library             *
  *         SCIP --- Solving Constraint Integer Programs                      *
  *                                                                           *
  *    Copyright (C) 2002-2021 Zuse Institute Berlin                          *
  *                                                                           *
  *                                                                           *
  *  SCIP is distributed under the terms of the ZIB Academic License.         *
  *                                                                           *
  *  You should have received a copy of the ZIB Academic License              *
  *  along with SCIP; see the file COPYING. If not visit scipopt.org.         *
  *                                                                           *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  @file   pricer_binpacking.c
  @brief  Binpacking variable pricer
  @author Timo Berthold
  @author Stefan Heinz

  The code is adapted to our terminology. Except for this, the code remains unchanged.
*/

/**@file   pricer_pattern.h
 * @brief  variable pricer to generate new patterns
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PRICER_PATTERN_H__
#define __PRICER_PATTERN_H__

#include "datapoints.h"

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the pattern variable pricer and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePricerPattern(
   SCIP*                 scip                /**< SCIP data structure */
   );

/** adds problem specific data to pricer and activates pricer */
SCIP_EXPORT
SCIP_RETCODE SCIPpricerPatternActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   absmax,             /**< maximum absolute value of coordinate in X */
   SCIP_CONS**           coverconss,         /**< covering constraints */
   int                   ncoverconss,        /**< number of covering constraints */
   SCIP_CONSHDLR*        conshdlr            /**< pointer to same/diff conshdlr */
   );


#ifdef __cplusplus
}
#endif

#endif
