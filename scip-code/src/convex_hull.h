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


/* constructs a CDD matrix containing points from a set of data points */
extern
dd_MatrixPtr constructGeneratorMatrixPoints(
   Datapoints*           datapoints,         /**< data points */
   int*                  selectedpoints,     /**< array of selected data points */
   int                   nselectedpoints     /**< number of selected data points */
   );

/** computes facet description of set of points */
extern
dd_MatrixPtr computeConvexHullFacets(
   SCIP*                 scip,               /**< SCIP instance */
   dd_MatrixPtr          generators,         /**< points generating convex hull */
   SCIP_Bool*            success             /**< pointer to store whether we were successful */
   );

#ifdef __cplusplus
}
#endif

#endif
