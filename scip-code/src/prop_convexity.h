/**@file   prop_convexity.h
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __PROP_CONVEXITY_H__
#define __PROP_CONVEXITY_H__


#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


/** creates the convexity propagator and includes it in SCIP */
SCIP_EXPORT
SCIP_RETCODE SCIPincludePropConvexity(
   SCIP*                 scip                /**< SCIP data structure */
   );

#ifdef __cplusplus
}
#endif

#endif
