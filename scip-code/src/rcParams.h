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
