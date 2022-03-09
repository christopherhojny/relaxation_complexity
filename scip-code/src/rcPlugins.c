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
