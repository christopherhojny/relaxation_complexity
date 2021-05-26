#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "math.h"

#include "datapoints.h"
#include "scip/scip.h"

#include "auxiliary_cdd.h"
#include "convex_hull.h"


dd_MatrixPtr constructGeneratorMatrixPoints(
   Datapoints*           datapoints,
   int*                  pointsincluster,
   int                   npointsincluster
   )
{
   dd_MatrixPtr generators;
   dd_colrange dim;
   dd_rowrange i;
   dd_colrange j;

   assert( datapoints != NULL );
   assert( pointsincluster != NULL || npointsincluster == 0 );
   assert( npointsincluster >= 0 );

   /* initialize generators matrix */
   dim = datapoints->dimension + 1;

   /* add points (we use homogeneous coordinates, i.e., first coordinate is 1) */
   if ( npointsincluster > 0 )
   {
      generators = dd_CreateMatrix(npointsincluster, dim);

      for (i = 0; i < npointsincluster; ++i)
      {
         dd_set_si(generators->matrix[i][0], 1);
         for (j = 1; j < dim; ++j)
            dd_set_si(generators->matrix[i][j], datapoints->points[pointsincluster[i]][j-1]);
      }
   }
   else
   {
      generators = dd_CreateMatrix(datapoints->ndatapoints, dim);

      for (i = 0; i < datapoints->ndatapoints; ++i)
      {
         dd_set_si(generators->matrix[i][0], 1);
         for (j = 1; j < dim; ++j)
            dd_set_si(generators->matrix[i][j], datapoints->points[i][j-1]);
      }
   }

   return generators;
}


dd_MatrixPtr computeConvexHullFacets(
   SCIP*                 scip,
   dd_MatrixPtr          generators,
   SCIP_Bool*            success
   )
{
   dd_PolyhedraPtr polyhedron;
   dd_MatrixPtr facets;
   dd_ErrorType err;

   assert( scip != NULL );
   assert( generators != NULL );
   assert( success != NULL );

   *success = FALSE;

   /* tell CDD that the generators are generators */
   generators->representation = dd_Generator;

   /* compute H-description */
   polyhedron = dd_DDMatrix2Poly(generators, &err);

   if ( err != dd_NoError )
   {
      SCIPdebugMsg(scip, "could not compute H-description\n");

      return NULL;
   }

   *success = TRUE;

   facets = dd_CopyInequalities(polyhedron);

   dd_FreePolyhedra(polyhedron);

   return facets;

}
