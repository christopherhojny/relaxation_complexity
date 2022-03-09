#ifndef __COMPUTE_SYMMETRIES_H__
#define __COMPUTE_SYMMETRIES_H__

#include "datapoints.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


extern
/** computes common symmetries of X and Y */
SCIP_RETCODE computeSymmetries(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int***                perms,              /**< pointer to store permutations */
   int*                  nperms,             /**< pointer to store number of permutations stored in perms */
   int*                  nmaxperms           /**< pointer to store maximum number of permutations fiiting in perms */
   );

#ifdef __cplusplus
}
#endif

#endif
