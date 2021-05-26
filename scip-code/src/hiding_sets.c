#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "math.h"

#include "convex_hull.h"
#include "datapoints.h"
#include "scip/scip.h"

#include "auxiliary_cdd.h"
#include "hiding_sets.h"


/** computes index sets of hiding set cuts, i.e., pairs of integer points forming a hiding set */
SCIP_RETCODE computeHidingSetCuts(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int**                 hidingsetidx,       /**< address to pointer for storing hiding set indices, i.e.,
                                              *   (tuples) of indices of infeasible points */
   int*                  nhidingsetidx,      /**< number of hiding set indices encoded in hidingsetidx */
   int*                  maxnhidingsetidx    /**< maximum number of entries in hidingsetidx */
   )
{
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   SCIP_Bool success;
   dd_rowrange rowidx;
   dd_colrange colidx;
   int yidx;
   int zidx;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( hidingsetidx != NULL );
   assert( nhidingsetidx != NULL );
   assert( maxnhidingsetidx != NULL );

   *nhidingsetidx = 0;
   *maxnhidingsetidx = 0;

   /* compute inequality description of conv(X) */
   dd_set_global_constants();
   generators = constructGeneratorMatrixPoints(X, NULL, 0);

   facetsconvexhull = computeConvexHullFacets(scip, generators, &success);

   dd_FreeMatrix(generators);

   if ( ! success )
   {
      dd_free_global_constants();
      return SCIP_OKAY;
   }

   /* for each pair of points (y,z) in Y, check whether they form a hiding set
    *
    * To this end, compute for each facet inequality a (lambda * y + (1-lambda)*z) <= b
    * <=> lambda a (y - z) <= b - a (z)
    * and solve it for lambda. Use these values to update the lower and upper bounds on
    * the convex multiplier. If the system is feasible, they form a hiding set.
    */

   /* iterate over (y,z) */
   for (yidx = 0; yidx < Y->ndatapoints; ++yidx)
   {
      for (zidx = yidx + 1; zidx < Y->ndatapoints; ++zidx)
      {
         SCIP_Real lblambda = 0.0;
         SCIP_Real ublambda = 1.0;

         /* iterate over facets */
         for (rowidx = 0; rowidx < facetsconvexhull->rowsize; ++rowidx)
         {
            SCIP_Real lhs = 0.0;
            SCIP_Real rhs;

            rhs = getReal(facetsconvexhull->matrix[rowidx][0]);

            for (colidx = 1; colidx < facetsconvexhull->colsize; ++colidx)
            {
               rhs += getReal(facetsconvexhull->matrix[rowidx][colidx]) * Y->points[zidx][colidx - 1];
               lhs -= getReal(facetsconvexhull->matrix[rowidx][colidx]) * Y->points[yidx][colidx - 1];
               lhs += getReal(facetsconvexhull->matrix[rowidx][colidx]) * Y->points[zidx][colidx - 1];
            }

            /* update lower bound or upper bound on lambda */
            if ( SCIPisLT(scip, lhs, 0.0) )
               lblambda = MAX(lblambda, rhs/lhs);
            else if ( SCIPisGT(scip, lhs, 0.0) )
               ublambda = MIN(ublambda, rhs/lhs);
            else if ( SCIPisLT(scip, rhs, 0.0) )
               break;

            if ( SCIPisGT(scip, lblambda, ublambda) )
               break;
         }

         /* check whether we have found a hiding set */
         if ( rowidx == facetsconvexhull->rowsize )
         {
            /* ensure that we can store all cuts */
            if ( *maxnhidingsetidx <= *nhidingsetidx + 2 )
            {
               if ( *maxnhidingsetidx == 0 )
               {
                  *maxnhidingsetidx = Y->ndatapoints;
                  SCIP_CALL( SCIPallocBlockMemoryArray(scip, hidingsetidx, *maxnhidingsetidx) );
               }
               else
               {
                  int newsize;

                  newsize = SCIPcalcMemGrowSize(scip, *maxnhidingsetidx + 2);
                  SCIP_CALL( SCIPreallocBlockMemoryArray(scip, hidingsetidx, *maxnhidingsetidx, newsize) );
                  *maxnhidingsetidx = newsize;
               }
            }

            (*hidingsetidx)[(*nhidingsetidx)++] = yidx;
            (*hidingsetidx)[(*nhidingsetidx)++] = zidx;
         }
      }
   }

   dd_FreeMatrix(facetsconvexhull);
   dd_free_global_constants();

   return SCIP_OKAY;
}
