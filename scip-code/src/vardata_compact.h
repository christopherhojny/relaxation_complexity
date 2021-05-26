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
