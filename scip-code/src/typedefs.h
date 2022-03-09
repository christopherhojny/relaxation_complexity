#ifndef __TYPEDEFS_H_
#define __TYPEDEFS_H_

#ifdef __cplusplus
extern "C" {
#endif

#define METHOD_COMPACT_MIP      0                   /**< compact MIP model */
#define METHOD_CG               1                   /**< column generation */
#define METHOD_HYBRID_COMPACT   2                   /**< column generation for generating bound and compact for rest */
#define METHOD_CONFLICT         3                   /**< conflict based model */
#define METHOD_HYBRID_CONFLICT  4                   /**< column generation for generating bound and conflict for rest */


#ifdef __cplusplus
}
#endif

#endif
