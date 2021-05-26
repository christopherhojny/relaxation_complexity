/**@file   rcPlugins.c
 * @brief  load SCIP plugins for computing RC
 * @author Christopher Hojny
 */

#include "rcPlugins.h"

#include "scip/scipdefplugins.h"


/** Include basic plugins needed for computing RC */
SCIP_RETCODE includeRCPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   SCIP_CALL( SCIPincludeDefaultPlugins(scip) );

   return SCIP_OKAY;
}
