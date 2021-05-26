/**@file   rcPlugins.h
 * @brief  load SCIP plugins for computing RC
 * @author Christopher Hojny
 */

#ifndef __RCPLUGINS_H__
#define __RCPLUGINS_H__

// SCIP include
#include <scip/scip.h>


#ifdef __cplusplus
extern "C" {
#endif

/** Include basic plugins needed for computing RC */
extern
SCIP_RETCODE includeRCPlugins(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
