/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program computeRC                             */
/*                                                                           */
/*    an implementation of a branch-and-cut and branch-and-price             */
/*    algorithm to compute the epsilon relaxation complexity of              */
/*    a full-dimensional lattice-convex set X and a finite set               */
/*    of points Y.                                                           */
/*                                                                           */
/*    Copyright (C) 2022-     Gennadiy Averkov, Christopher Hojny,           */
/*                            Matthias Schymura                              */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/*    Copyright (C) 2002-2022 Zuse Institute Berlin                          */
/*                                                                           */
/*       mailto: scip@zib.de                                                 */
/*       Licensed under the Apache License, Version 2.0                      */
/*                                                                           */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   convex_hull.c
 * @brief  file implementing functions for convex hull computations
 * @author Christopher Hojny
 */

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "math.h"

#include "datapoints.h"
#include "scip/scip.h"

#include "auxiliary_cdd.h"
#include "convex_hull.h"


/** constructs a CDD matrix containing points from a set of data points */
dd_MatrixPtr constructGeneratorMatrixPoints(
   Datapoints*           datapoints,         /**< data points */
   int*                  selectedpoints,     /**< array of selected data points */
   int                   nselectedpoints     /**< number of selected data points */
   )
{
   dd_MatrixPtr generators;
   dd_colrange dim;
   dd_rowrange i;
   dd_colrange j;

   assert( datapoints != NULL );
   assert( selectedpoints != NULL || nselectedpoints == 0 );
   assert( nselectedpoints >= 0 );

   /* initialize generators matrix */
   dim = datapoints->dimension + 1;

   /* add points (we use homogeneous coordinates, i.e., first coordinate is 1) */
   if ( nselectedpoints > 0 )
   {
      generators = dd_CreateMatrix(nselectedpoints, dim);

      for (i = 0; i < nselectedpoints; ++i)
      {
         dd_set_si(generators->matrix[i][0], 1);
         for (j = 1; j < dim; ++j)
            dd_set_si(generators->matrix[i][j], datapoints->points[selectedpoints[i]][j-1]);
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


/** computes facet description of set of points */
dd_MatrixPtr computeConvexHullFacets(
   SCIP*                 scip,               /**< SCIP instance */
   dd_MatrixPtr          generators,         /**< points generating convex hull */
   SCIP_Bool*            success             /**< pointer to store whether we were successful */
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

   facets = dd_CopyInequalities(polyhedron);

   if ( set_card(facets->linset) > 0 )
   {
      SCIPdebugMsg(scip, "polyhedron is not full-dimensional");

      return NULL;
   }

   *success = TRUE;

   dd_FreePolyhedra(polyhedron);

   return facets;

}
