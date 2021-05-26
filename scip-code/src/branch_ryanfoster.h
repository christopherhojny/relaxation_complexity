/**@file   branch_ryanfoster.h
 * @brief  Ryan/Foster branching rule
 * @author Christopher Hojny
 */

#ifndef __SCIP_BRANCH_RYANFOSTER_H__
#define __SCIP_BRANCH_RYANFOSTER_H__


#include "scip/scip.h"

/** creates the ryanfoster branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanFoster(
   SCIP*                 scip                /**< SCIP data structure */
   );

#endif
