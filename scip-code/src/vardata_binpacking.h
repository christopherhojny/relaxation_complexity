/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*    Copyright (C) 2002-2020 Zuse Institute Berlin                          */
/*                                                                           */
/*                                                                           */
/*  SCIP is distributed under the terms of the ZIB Academic License.         */
/*                                                                           */
/*  You should have received a copy of the ZIB Academic License              */
/*  along with SCIP; see the file COPYING. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   vardata_binpacking.h
 * @brief  Variable data containing the ids of constraints in which the variable appears
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Christopher Hojny
 *
 * This file implements the handling of the variable data which is attached to each file. See SCIP_VarData and \ref BINPACKING_PRICER.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_VARDATA_BINPACKING__
#define __SCIP_VARDATA_BINPACKING__

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** create variable data */
SCIP_RETCODE SCIPvardataCreateBinpacking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int*                  consids,            /**< array of constraints ids */
   int                   nconss,             /**< number of constraints */
   SCIP_Real*            inequality,         /**< array encoding inequality defining the set */
   int                   nentries            /**< number of entries in inequality */
   );

/** get number of constraints */
int SCIPvardataGetNConsids(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns sorted constraint id array */
int* SCIPvardataGetConsids(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** get number of constraints */
int SCIPvardataGetNEntries(
   SCIP_VARDATA*         vardata             /**< variable data */
   );

/** returns sorted constraint id array */
SCIP_Real* SCIPvardataGetInequality(
   SCIP_VARDATA*         vardata             /**< variable data */
   );


/** creates variable */
SCIP_RETCODE SCIPcreateVarBinpacking(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_Bool             initial,            /**< should var's column be present in the initial root LP? */
   SCIP_Bool             removable,          /**< is var's column removable from the LP (due to aging or cleanup)? */
   SCIP_Bool             deletable,          /**< is var's column deletable from the LP (due to aging or cleanup)? */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   );

#ifdef __cplusplus
}
#endif

#endif
