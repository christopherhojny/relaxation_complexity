#ifndef __CONVEX_HULL_H__
#define __CONVEX_HULL_H__

#include "datapoints.h"
#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif


extern
dd_MatrixPtr constructGeneratorMatrixPoints(
   Datapoints*           datapoints,
   int*                  pointsincluster,
   int                   npointsincluster
   );


extern
dd_MatrixPtr computeConvexHullFacets(
   SCIP*                 scip,
   dd_MatrixPtr          generators,
   SCIP_Bool*            success
   );

#ifdef __cplusplus
}
#endif

#endif
