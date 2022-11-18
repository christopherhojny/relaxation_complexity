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
 * @brief  implements the setup of the basic models
 * @author Christopher Hojny
 */

#include <stdio.h>
#include <time.h>

#include <scip/scip.h>
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"

#include "rcParams.h"
#include "rcPlugins.h"

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"
#include "math.h"
#include "convex_hull.h"

#include "auxiliary_cdd.h"
#include "branch_ryanfoster.h"
#include "cons_conflict.h"
#include "cons_samediff.h"
#include "maximal_separation.h"
#include "pricer_pattern.h"
#include "problem_rc.h"
#include "probdata_rc_cg.h"
#include "probdata_rc_compact.h"
#include "probdata_rc_conflict.h"
#include "prop_convexity.h"
#include "typedefs.h"
#include "vardata_binpacking.h"


/** removes points from the set of infeasible points that are no observers */
static
SCIP_RETCODE SCIPremoveNonobservers(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y                   /**< set of infeasible points */
   )
{
   Datapoints* testset;
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   SCIP_Real val;
   int p;
   int q;
   int i;
   int f;
   int npoints;
   int npointsold;
   int nremoved = 0;
   SCIP_Bool success = TRUE;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );

   /* create datapoints containing points in X and one additional point from Y */
   SCIP_CALL( SCIPallocBlockMemory(scip, &testset) );
   testset->ndatapoints = X->ndatapoints + 1;
   testset->dimension = X->dimension;

   npointsold = Y->ndatapoints;

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &testset->points, testset->ndatapoints) );
   for (p = 0; p < testset->ndatapoints; ++p)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &testset->points[p], testset->dimension) );
   }

   /* copy X */
   for (p = 0; p < X->ndatapoints; ++p)
   {
      for (i = 0; i < testset->dimension; ++i)
         testset->points[p][i] = X->points[p][i];
   }

   /* for each point in Y, check whether y is an observer w.r.t. (X,Y) */
   npoints = Y->ndatapoints - 1;
   while ( npoints >= 0 )
   {
      p = npoints;

      for (i = 0; i < Y->dimension; ++i)
         testset->points[X->ndatapoints][i] = Y->points[p][i];

      /* compute convex hull of X and y */
      dd_set_global_constants();

      generators = constructGeneratorMatrixPoints(testset, NULL, 0);
      facetsconvexhull = computeConvexHullFacets(scip, generators, &success);

      if ( ! success )
      {
         dd_free_global_constants();

         return SCIP_ERROR;
      }

      dd_FreeMatrix(generators);

      /* check whether the convex hull contains other points from Y */
      for (q = 0; q < Y->ndatapoints; ++q)
      {
         /* skip the point itself */
         if ( q == p )
            continue;

         /* iterate over facets and check whether q satisfies all of them */
         for (f = 0; f < facetsconvexhull->rowsize; ++f)
         {
            val = getReal(facetsconvexhull->matrix[f][0]);
            for (i = 1; i <= Y->dimension; ++i)
               val += getReal(facetsconvexhull->matrix[f][i]) * Y->points[q][i - 1];

            /* inequality is strictly violated */
            if ( SCIPisLT(scip, val, 0.0) )
               break;
         }

         /* point q is contained in the convex hull: STOP,
            remove p from list, because it cannot be an observer
         */
         if ( f == facetsconvexhull->rowsize )
            break;
      }

      /* possible remove p */
      if ( q != Y->ndatapoints )
      {
         ++nremoved;

         for (i = npoints; i < Y->ndatapoints - 1; ++i)
         {
            SCIP_Real* point;
            point = Y->points[i];
            Y->points[i] = Y->points[i+1];
            Y->points[i+1] = point;
         }
         Y->ndatapoints -= 1;
      }

      dd_FreeMatrix(facetsconvexhull);
      dd_free_global_constants();

      --npoints;
   }

   /* possibly reallocate memory in Y */
   if ( npointsold > Y->ndatapoints )
   {
      for (p = npointsold - 1; p >= Y->ndatapoints; --p)
      {
         SCIPfreeBlockMemoryArray(scip, &Y->points[p], Y->dimension);
      }
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &Y->points, npointsold, Y->ndatapoints) );
   }

   /* free memory */
   for (p = X->ndatapoints; p >= 0; --p)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &testset->points[p], testset->dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &testset->points, testset->ndatapoints);
   SCIPfreeBlockMemory(scip, &testset);

   SCIPinfoMessage(scip, NULL, "PRESOLVE: removed %d points from the set of infeasible points that are not observers\n\n",
      nremoved);

   return SCIP_OKAY;
}


/** writes a solution of the CG model to a file */
static
SCIP_RETCODE writeSolutionCG(
   SCIP*                 scip,               /**< SCIP instance */
   char*                 filename            /**< file name */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   FILE* fp;
   int i;
   int j;
   SCIP_Bool display = FALSE;

   assert( scip != NULL );

   /* get best solution */
   sol = SCIPgetBestSol(scip);

   if ( sol == NULL )
   {
      printf("No solution could be found\n");
      return SCIP_OKAY;
   }

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   if ( filename == 0 )
   {
      printf("\nBest solution found:\n");
      display = TRUE;
   }
   else
      fp = fopen(filename, "w");

   /* iterate over variables and check which sets have been selected in solution */
   for (i = 0; i < SCIPprobdataGetNVars(probdata); ++i)
   {
      /* return inequality defining used set */
      if ( SCIPgetSolVal(scip, sol, SCIPprobdataGetVars(probdata)[i]) > 0.5 )
      {
         SCIP_VAR* var;
         SCIP_VARDATA* vardata;

         var = SCIPprobdataGetVars(probdata)[i];
         vardata = SCIPvarGetData(var);
         display ? printf("%f", SCIPvardataGetInequality(vardata)[0]) : fprintf(fp, "%f", SCIPvardataGetInequality(vardata)[0]);
         for (j = 1; j < SCIPvardataGetNEntries(vardata); ++j)
         {
            display ? printf(" %f", SCIPvardataGetInequality(vardata)[j]) : fprintf(fp, " %f", SCIPvardataGetInequality(vardata)[j]);
         }
         display ? printf("\n") : fprintf(fp, "\n");
      }
   }

   if ( ! display )
      fclose(fp);

   return SCIP_OKAY;
}


/** writes a solution of the compact model to a file */
static
SCIP_RETCODE writeSolutionCompact(
   SCIP*                 scip,               /**< SCIP instance */
   char*                 filename,           /**< file name */
   int                   ninequalities,      /**< number of possible inequalities in formulation */
   int                   dimension           /**< dimension of data points */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   SCIP_VAR*** lhsvars;
   SCIP_VAR** rhsvars;
   SCIP_VAR* var;
   FILE* fp;
   int i;
   int j;
   SCIP_Bool display = FALSE;

   assert( scip != NULL );

   /* get best solution */
   sol = SCIPgetBestSol(scip);

   if ( sol == NULL )
   {
      printf("No solution could be found\n");
      return SCIP_OKAY;
   }

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   if ( filename == 0 )
   {
      display = TRUE;
      printf("\nBest solution found:\n");
   }
   else
      fp = fopen(filename, "w");

   lhsvars = SCIPprobdataGetLhsvars(probdata);
   rhsvars = SCIPprobdataGetRhsvars(probdata);

   /* iterate over variables and check which sets have been selected in solution */
   for (i = 0; i < ninequalities; ++i)
   {
      /* return inequality defining used set */
      if ( SCIPgetSolVal(scip, sol, SCIPprobdataGetUsedvars(probdata)[i]) > 0.5 )
      {
         var = rhsvars[i];
         display ? printf("%f", SCIPgetSolVal(scip, sol, var)) : fprintf(fp, "%f", SCIPgetSolVal(scip, sol, var));

         for (j = 0; j < dimension; ++j)
         {
            var = lhsvars[i][j];
            display ? printf(" %f", -SCIPgetSolVal(scip, sol, var)) : fprintf(fp, " %f", -SCIPgetSolVal(scip, sol,var));
         }
         display ? printf("\n") : fprintf(fp, "\n");
      }
   }

   if ( ! display )
      fclose(fp);

   return SCIP_OKAY;
}


/** computes inequality separating a set of points */
static
SCIP_RETCODE SCIPcomputeInequality(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  selectedpoints,     /**< array of selected points from Y */
   int                   nselectedpoints,    /**< number of points stored in selectedpoints */
   SCIP_Real*            inequality,         /**< allocated array to store inequality of type b + ax >=0 */
   int                   dim,                /**< dimension of points in X and Y */
   int                   absmaxX             /**< maximum absolute coordinate of a point in X */
   )
{
   SCIP* subscip;
   SCIP_VAR** ineqvars;
   SCIP_CONS** validconss;
   SCIP_CONS** sepaconss;
   SCIP_SOL* sol;
   SCIP_Real eps;
   int i;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( selectedpoints != NULL );
   assert( nselectedpoints > 0 );
   assert( inequality != NULL );
   assert( dim > 0 );

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );
   SCIP_CALL( SCIPallocBufferArray(scip, &ineqvars, dim + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &validconss, X->ndatapoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaconss, nselectedpoints) );

   SCIP_CALL( SCIPcreate(&subscip) );

   SCIP_CALL( SCIPcreateSeparationLP(subscip, X, Y, selectedpoints, nselectedpoints, eps, absmaxX,
         ineqvars, validconss, sepaconss) );

   SCIP_CALL( SCIPsolve(subscip) );

   sol = SCIPgetBestSol(subscip);
   assert( sol != NULL );

   inequality[dim] = SCIPgetSolVal(subscip, sol, ineqvars[dim]);
   for (i = 0; i < dim; ++i)
      inequality[i] = - SCIPgetSolVal(subscip, sol, ineqvars[i]);

   SCIP_CALL( SCIPfreeSeparationLP(subscip, X->ndatapoints, nselectedpoints, dim, ineqvars, validconss, sepaconss) );
   SCIP_CALL( SCIPfree(&subscip) );
   SCIPfreeBufferArray(scip, &sepaconss);
   SCIPfreeBufferArray(scip, &validconss);
   SCIPfreeBufferArray(scip, &ineqvars);

   return SCIP_OKAY;
}


/** writes a solution of the conflict model to a file */
static
SCIP_RETCODE writeSolutionConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   char*                 filename,           /**< file name */
   int                   ub,                 /**< upper bound on RC(X,Y) */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   absmaxX             /**< maximum absolute coordinate of a point in X */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_SOL* sol;
   FILE* fp;
   SCIP_VAR** isuedvars;
   SCIP_VAR*** violatedvars;
   SCIP_Real* inequality;
   int* violatedpoints;
   int nviolatedpoints;
   int nY;
   int dim;
   int i;
   int j;
   SCIP_Bool display = FALSE;

   assert( scip != NULL );
   assert( ub > 0 );
   assert( X != NULL );
   assert( Y != NULL );

   /* get best solution */
   sol = SCIPgetBestSol(scip);

   if ( sol == NULL )
   {
      printf("No solution could be found\n");
      return SCIP_OKAY;
   }

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   isuedvars = SCIPprobdataConflictGetIsusedvars(probdata);
   violatedvars = SCIPprobdataConflictGetViolatedvars(probdata);
   nY = Y->ndatapoints;
   dim = Y->dimension;

   assert( isuedvars != NULL );
   assert( violatedvars != NULL );

   if ( filename == 0 )
   {
      display = TRUE;
      printf("\nBest solution found:\n");
   }
   else
      fp = fopen(filename, "w");

   SCIP_CALL( SCIPallocBufferArray(scip, &violatedpoints, nY) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inequality, dim + 1) );

   /* for each used inequality, compute the inequality and return it */
   for (i = 0; i < ub; ++i)
   {
      /* skip unused inequalities */
      if ( SCIPgetSolVal(scip, sol, isuedvars[i]) < 0.5 )
         continue;

      /* get points separated by inequality */
      nviolatedpoints = 0;
      for (j = 0; j < nY; ++j)
      {
         if ( SCIPgetSolVal(scip, sol, violatedvars[j][i]) > 0.5 )
            violatedpoints[nviolatedpoints++] = j;
      }

      SCIP_CALL( SCIPcomputeInequality(scip, X, Y, violatedpoints, nviolatedpoints, inequality, dim, absmaxX) );

      display ? printf("%f", inequality[dim]) : fprintf(fp, "%f", inequality[dim]);
      for (j = 0; j < dim; ++j)
      {
         display ? printf(" %f", inequality[j]) : fprintf(fp, " %f", inequality[j]);
      }
      display ? printf("\n") : fprintf(fp, "\n");
   }

   SCIPfreeBufferArray(scip, &inequality);
   SCIPfreeBufferArray(scip, &violatedpoints);

   return SCIP_OKAY;
}


/** creates initial model for computing RC */
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX,            /**< maximum absolute value of a coordinate in X */
   SCIP_Bool             secondphase,        /**< whether we are in the second phase of the hybrid approach */
   SCIP_Real*            heurtime,           /**< pointer to store time spent in heuristics */
   SCIP_Bool*            success             /**< pointer to store whether model could be built successfully */
   )
{
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   SCIP_Real timelimit;
   int method;
   int allocub;
   int i;
   int j;
   int k;
   SCIP_Bool dopreprocessing;
   SCIP_Real** inequalities = NULL;
   int** separatedpoints = NULL;
   int* nseparatedpoints = NULL;
   int ninequalities = -1;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub != NULL );
   assert( *ub > 0 );
   assert( heurtime != NULL || secondphase );
   assert( success != NULL );

   *success = TRUE;

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/preprocessing", &dopreprocessing) );
   SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );

   if ( ! secondphase && dopreprocessing )
   {
      SCIP_CALL( SCIPremoveNonobservers(scip, X, Y) );
   }

   SCIP_CALL( SCIPgetIntParam(scip, "rc/method", &method) );

   if ( ! secondphase )
   {
      clock_t starttime;
      clock_t endtime;
      SCIP_Real usedtime;
      SCIP_Bool ubcomputed;
      int tmpninequalities;
      SCIP_Real eps;

      SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

      /*
       * compute initial solutions
       */
      *heurtime = 0.0;

      SCIPinfoMessage(scip, NULL, "   compute facet description of conv(X)\n");
      starttime = clock();

      /* facet description */
      dd_set_global_constants();
      generators = constructGeneratorMatrixPoints(X, NULL, 0);

      facetsconvexhull = computeConvexHullFacets(scip, generators, &ubcomputed);

      dd_FreeMatrix(generators);
      endtime = clock();
      usedtime = (SCIP_Real) (endtime - starttime) / CLOCKS_PER_SEC;
      *heurtime += usedtime;
      timelimit -= usedtime;

      if ( SCIPisLT(scip, timelimit, 0.0) )
      {
         SCIPinfoMessage(scip, NULL, "time limit hit during computation of facets of conv(X)\n");
         *success = FALSE;

         return SCIP_OKAY;
      }

      if ( ubcomputed && facetsconvexhull->rowsize < *ub)
      {
         *ub = facetsconvexhull->rowsize;
         SCIPinfoMessage(scip, NULL, "   facet computation successful, found heuristic solution with %d inequalities in %.1fs\n",
            *ub, usedtime);
      }
      else
      {
         dd_free_global_constants();
         ubcomputed = FALSE;
         SCIPinfoMessage(scip, NULL, "   facet computation finished, did not find better solution in %.1fs\n",
            usedtime);
      }

      /* prepare memory to store greedy solution */
      *ub = MIN(*ub, Y->ndatapoints);
      allocub = *ub;

      SCIP_CALL( SCIPallocBufferArray(scip, &inequalities, *ub) );
      for (i = 0; i < *ub; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &inequalities[i], X->dimension + 1) );
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &separatedpoints, *ub) );
      for (i = 0; i < *ub; ++i)
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &separatedpoints[i], Y->ndatapoints) );
      }
      SCIP_CALL( SCIPallocBufferArray(scip, &nseparatedpoints, *ub) );

#if 0
      /* compute greedy solution based on standard order */
      SCIPinfoMessage(scip, NULL, "   compute greedy solution based on standard order\n");
      starttime = clock();

      ninequalities = -1;
      tmpninequalities = -1;
      SCIP_CALL( heurGreedy(scip, X, Y, X->ndatapoints, Y->ndatapoints, X->dimension, absmaxX, *ub,
            GREEDY_SORT_STANDARD, method, inequalities, separatedpoints, nseparatedpoints, &tmpninequalities,
            ubcomputed ? facetsconvexhull : NULL, timelimit) );
      endtime = clock();
      usedtime = (SCIP_Real) (endtime - starttime) / CLOCKS_PER_SEC;
      *heurtime += usedtime;
      timelimit -= usedtime;

      if ( tmpninequalities > 0 )
      {
         /* we have found a better solution, because tmpninequalities is positive */
         SCIPinfoMessage(scip, NULL, "   standard greedy found solution with %d inequalities in %.1fs\n",
            tmpninequalities, usedtime);
         *ub = tmpninequalities;
         ninequalities = tmpninequalities;
      }
      else
         SCIPinfoMessage(scip, NULL, "   standard greedy did not find a better solution in %.1fs\n",
            usedtime);

      if ( SCIPisLT(scip, timelimit, 0.0) )
      {
         SCIPinfoMessage(scip, NULL, "time limit hit during standard greedy\n");
         *success = FALSE;

         goto FREEMEMORY;
      }
#endif

      SCIPinfoMessage(scip, NULL, "   compute greedy solution based on lattice distance order\n");
      starttime = clock();

      tmpninequalities = -1;
      SCIP_CALL( heurGreedy(scip, X, Y, X->ndatapoints, Y->ndatapoints, X->dimension, absmaxX, *ub,
            GREEDY_SORT_DISTANCE, method, inequalities, separatedpoints, nseparatedpoints, &tmpninequalities,
            ubcomputed ? facetsconvexhull : NULL, timelimit) );
      endtime = clock();
      usedtime = (SCIP_Real) (endtime - starttime) / CLOCKS_PER_SEC;
      *heurtime += usedtime;
      timelimit -= usedtime;

      if ( tmpninequalities > 0 )
      {
         /* we have found a better solutoin, because tmpninequalities is positive */
         SCIPinfoMessage(scip, NULL, "   lattice distance greedy found solution with %d inequalities in %.1fs\n",
            tmpninequalities, usedtime);
         *ub = tmpninequalities;
         ninequalities = tmpninequalities;
      }
      else
         SCIPinfoMessage(scip, NULL, "   lattice distance greedy did not find a better solution in %.1fs\n",
            usedtime);

      if ( SCIPisLT(scip, timelimit, 0.0) )
      {
         SCIPinfoMessage(scip, NULL, "time limit hit during standard greedy\n");
         *success = FALSE;

         goto FREEMEMORY;
      }

      /* store facet solution if it is better than the best greedy solution */
      if ( ubcomputed && facetsconvexhull->rowsize <= *ub
         && (facetsconvexhull->rowsize < ninequalities || ninequalities < 0) )
      {
         for (i = 0; i < facetsconvexhull->rowsize; ++i)
         {
            for (j = 0; j < facetsconvexhull->colsize; ++j)
               inequalities[i][j] = getReal(facetsconvexhull->matrix[i][j]);
         }
         ninequalities = facetsconvexhull->rowsize;

         for (i = 0; i < ninequalities; ++i)
         {
            nseparatedpoints[i] = 0;
            for (j = 0; j < Y->ndatapoints; ++j)
            {
               SCIP_Real viol;

               viol = inequalities[i][0];
               for (k = 0; k < X->dimension; ++k)
                  viol += Y->points[j][k] * inequalities[i][k + 1];

               if ( SCIPisLE(scip, viol, -eps) )
                  separatedpoints[i][nseparatedpoints[i]++] = j;
            }
         }
      }
   }

   if ( method == METHOD_COMPACT_MIP || (method == METHOD_HYBRID_COMPACT && secondphase) )
   {
      SCIP_CALL( SCIPprobdataCreateCompact(scip, "name", X, Y, ub, lb, absmaxX,
            inequalities, separatedpoints, nseparatedpoints, ninequalities, allocub) );
   }
   else if ( method == METHOD_CONFLICT || (method == METHOD_HYBRID_CONFLICT && secondphase) )
   {
      /* SCIP_CALL( SCIPprobdataCreateConflict(scip, "name", X, Y, ub, lb, absmaxX) ); */
      SCIP_CALL( SCIPprobdataCreateConflict(scip, "name", X, Y, ub, lb, absmaxX,
            inequalities, separatedpoints, nseparatedpoints, ninequalities, allocub) );
   }
   else
   {
      SCIP_CALL( SCIPincludePricerPattern(scip) );
      SCIP_CALL( SCIPincludeBranchruleRyanFoster(scip) );
      SCIP_CALL( SCIPincludeConshdlrSamediff(scip) );
      SCIP_CALL( SCIPsetBoolParam(scip, "pricing/delvars", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "propagating/convexity/enabled", FALSE) );

      /* we have to ensure that SCIP does not remove the singleton sets, otherwise, the problem
       * might become infeasible after branching
       */
      if ( method == METHOD_CG )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "constraints/setppc/maxprerounds", 0) );
      }

      SCIP_CALL( SCIPprobdataCreateCG(scip, "name", X, Y, absmaxX,
            inequalities, separatedpoints, nseparatedpoints, ninequalities, allocub) );
   }

   if ( ! secondphase )
   {
   FREEMEMORY:
      SCIPfreeBufferArray(scip, &nseparatedpoints);
      for (i = allocub - 1; i >= 0; --i)
      {
         SCIPfreeBufferArray(scip, &separatedpoints[i]);
      }
      SCIPfreeBufferArray(scip, &separatedpoints);
      for (i = allocub - 1; i >= 0; --i)
      {
         SCIPfreeBufferArray(scip, &inequalities[i]);
      }
      SCIPfreeBufferArray(scip, &inequalities);
   }

   return SCIP_OKAY;
}

/** free problem data */
SCIP_RETCODE SCIPfreeModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y                   /**< pointer to data points of Y */
   )
{
   int i;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );

   /* free points encoded in Y */
   for (i = Y->ndatapoints - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &Y->points[i], Y->dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &Y->points, Y->ndatapoints);
   SCIPfreeBlockMemory(scip, &Y);

   /* free points encoded in X */
   for (i = X->ndatapoints - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &X->points[i], X->dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &X->points, X->ndatapoints);
   SCIPfreeBlockMemory(scip, &X);

   return SCIP_OKAY;
}

/** solves the problem */
SCIP_RETCODE SCIPsolveRCproblem(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   ub,                 /**< upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX,            /**< maximum absolute value of a coordinate in X */
   char*                 filename,           /**< name of file for storing solution (or NULL) */
   SCIP_Real             timelimit           /**< time limit */
   )
{
   int method;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub > 0 );

   SCIP_CALL( SCIPgetIntParam(scip, "rc/method", &method) );

   if ( method == METHOD_HYBRID_COMPACT )
   {
      SCIP* newscip;
      SCIP_PROBDATA* newprobdata;
      SCIP_PROBDATA* oldprobdata;
      SCIP_SOL* newsol;
      SCIP_SOL* oldsol;
      SCIP_VAR* solvar;
      SCIP_VARDATA* vardata;
      SCIP_Real dualbound;
      SCIP_Real primalbound;
      int* usedinequalities;
      int nusedinequalities;
      int newlb;
      int newub;
      int i;
      int j;
      int k;
      SCIP_Bool stored;
      SCIP_Bool success;

      SCIP_VAR*** lhsvars;
      SCIP_VAR** rhsvars;
      SCIP_VAR*** violatedvars;
      SCIP_VAR** isusedvars;
      SCIP_CONS*** indicatorconss;
      SCIP_Real* inequality;
      SCIP_VAR* var;
      SCIP_Real slackval;
      SCIP_Real newtimelimit;
      SCIP_Real eps;
      int* pattern;
      int dimension;
      int len;

      /* solve the root node relaxation */
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
      SCIP_CALL( SCIPsetBoolParam(scip, "propagating/convexity/enabled", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "propagating/intersection/enabled", FALSE) );
      SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts", 0) );
      SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

      SCIP_CALL( SCIPsolve(scip) );

      /* output statistics */
      SCIPinfoMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );

      dualbound = SCIPgetDualbound(scip);
      primalbound = SCIPgetPrimalbound(scip);

      if ( SCIPisEQ(scip, primalbound, MAX(lb, dualbound)) )
      {
         SCIP_CALL( writeSolutionCG(scip, filename) );
         return SCIP_OKAY;
      }

      newtimelimit = timelimit - SCIPgetSolvingTime(scip);
      if ( SCIPisLE(scip, newtimelimit, 0.0) )
         return SCIP_OKAY;

      /* create new SCIP instance for compact model */
      SCIP_CALL( SCIPcreate(&newscip) );

      /* load basic plugins*/
      SCIP_CALL( includeRCPlugins(newscip) );

      /* add our own parameters */
      SCIP_CALL( setSCIPParameters(newscip) );
      SCIP_CALL( addRCParameters(newscip) );
      SCIP_CALL( SCIPcopyParamSettings(scip, newscip) );
      SCIP_CALL( SCIPsetIntParam(newscip,"presolving/maxrestarts", -1) );
      SCIP_CALL( SCIPsetRealParam(newscip, "limits/time", newtimelimit) );
      SCIP_CALL( SCIPsetLongintParam(newscip, "limits/nodes", -1) );

      newlb = (int) MAX(SCIPceil(scip, dualbound), lb);
      newub = (int) MIN(SCIPfloor(scip, primalbound), ub);

      SCIP_CALL( SCIPcreateModel(newscip, X, Y, &newub, newlb, absmaxX, TRUE, NULL, &success) );
      assert( success );

      /* copy solution */
      newprobdata = SCIPgetProbData(newscip);
      assert( newprobdata != NULL );

      oldprobdata = SCIPgetProbData(scip);
      assert( oldprobdata != NULL );

      oldsol = SCIPgetBestSol(scip);
      SCIP_CALL( SCIPcreateOrigSol(newscip, &newsol, NULL) );

      violatedvars = SCIPprobdataGetViolatedvars(newprobdata);
      isusedvars = SCIPprobdataGetUsedvars(newprobdata);
      lhsvars = SCIPprobdataGetLhsvars(newprobdata);
      rhsvars = SCIPprobdataGetRhsvars(newprobdata);
      indicatorconss = SCIPprobdataGetLinkviolconss(newprobdata);

      /* copy solution values */
      if ( SCIPisGE(scip, newub, SCIPround(scip, primalbound)) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &usedinequalities, (int) SCIPround(scip, primalbound)) );
         nusedinequalities = 0;
         dimension = X->dimension;

         SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

         for (i = 0; i < SCIPprobdataGetNVars(oldprobdata); ++i)
         {
            solvar = SCIPprobdataGetVars(oldprobdata)[i];
            assert( solvar != NULL );

            /* skip unused variables */
            if ( SCIPgetSolVal(scip, oldsol, solvar) < 0.5 )
               continue;

            usedinequalities[nusedinequalities] = i;
            vardata = SCIPvarGetData(solvar);
            assert( vardata != NULL );

            /* copy inequality */
            inequality = SCIPvardataGetInequality(vardata);
            SCIP_CALL( SCIPsetSolVal(newscip, newsol, rhsvars[nusedinequalities], inequality[0]) );
            for (j = 0; j < dimension; ++j)
            {
               SCIP_CALL( SCIPsetSolVal(newscip, newsol, lhsvars[nusedinequalities][j], -inequality[j+1]) );
            }

            /* set that inequality is used */
            SCIP_CALL( SCIPsetSolVal(newscip, newsol, isusedvars[nusedinequalities], 1.0) );

            /* set cut pattern */
            pattern = SCIPvardataGetConsids(vardata);
            len = SCIPvardataGetNConsids(vardata);
            for (j = 0; j < len; ++j)
            {
               SCIP_CALL( SCIPsetSolVal(newscip, newsol, violatedvars[pattern[j]][nusedinequalities], 1.0) );
            }
            ++nusedinequalities;
         }
         assert( nusedinequalities == SCIPround(scip, primalbound) );

         /* set values of slack variables in indicator constraints */
         for (i = 0; i < Y->ndatapoints; ++i)
         {
            for (j = 0; j < nusedinequalities; ++j)
            {
               var = SCIPgetSlackVarIndicator(indicatorconss[i][j]);
               slackval = eps;

               solvar = SCIPprobdataGetVars(oldprobdata)[usedinequalities[j]];
               vardata = SCIPvarGetData(solvar);
               inequality = SCIPvardataGetInequality(vardata);

               /* compute -ay + b: note, stored inequality is already in shape b - ax >= 0 */
               for (k = 0; k < dimension; ++k)
               {
                  slackval += inequality[k+1] * Y->points[i][k];
               }
               slackval += inequality[0];

               SCIP_CALL( SCIPsetSolVal(newscip, newsol, var, MAX(slackval, 0.0)) );
            }
         }

         SCIPfreeBufferArray(scip, &usedinequalities);

         SCIP_CALL( SCIPaddSol(newscip, newsol, &stored) );
         SCIP_CALL( SCIPfreeSol(newscip, &newsol) );
      }

      SCIP_CALL( SCIPsetBoolParam(newscip, "propagating/convexity/enabled", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(newscip, "propagating/intersection/enabled", TRUE) );
      SCIP_CALL( SCIPsolve(newscip) );
      SCIP_CALL( SCIPprintStatistics(newscip, NULL) );

      SCIP_CALL( writeSolutionCompact(newscip, filename, SCIPround(newscip, SCIPgetPrimalbound(newscip)), X->dimension) );

      /* free new structures */
      SCIP_CALL( SCIPfreeTransform(newscip) );
      SCIP_CALL( SCIPfree(&newscip) );
   }
   else if ( method == METHOD_HYBRID_CONFLICT )
   {
      SCIP* newscip;
      SCIP_PROBDATA* newprobdata;
      SCIP_PROBDATA* oldprobdata;
      SCIP_SOL* newsol;
      SCIP_SOL* oldsol;
      SCIP_VAR* solvar;
      SCIP_VARDATA* vardata;
      SCIP_Real dualbound;
      SCIP_Real primalbound;
      int* usedinequalities;
      int nusedinequalities;
      int newlb;
      int newub;
      int i;
      int j;
      SCIP_Bool stored;
      SCIP_Bool success;

      SCIP_VAR*** violatedvars;
      SCIP_VAR** isusedvars;
      SCIP_Real newtimelimit;
      int* pattern;
      int len;

      /* solve the root node relaxation */
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
      SCIP_CALL( SCIPsetBoolParam(scip, "propagating/convexity/enabled", FALSE) );
      SCIP_CALL( SCIPsetBoolParam(scip, "propagating/intersection/enabled", FALSE) );
      SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts", 0) );
      SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );

      SCIP_CALL( SCIPsolve(scip) );

      /* output statistics */
      SCIPinfoMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );

      dualbound = SCIPgetDualbound(scip);
      primalbound = SCIPgetPrimalbound(scip);

      if ( SCIPisEQ(scip, primalbound, MAX(lb, dualbound)) )
      {
         SCIP_CALL( writeSolutionCG(scip, filename) );
         return SCIP_OKAY;
      }

      newtimelimit = timelimit - SCIPgetSolvingTime(scip);
      if ( SCIPisLE(scip, newtimelimit, 0.0) )
         return SCIP_OKAY;

      /* create new SCIP instance for compact model */
      SCIP_CALL( SCIPcreate(&newscip) );

      /* load basic plugins*/
      SCIP_CALL( includeRCPlugins(newscip) );

      /* add our own parameters */
      SCIP_CALL( setSCIPParameters(newscip) );
      SCIP_CALL( addRCParameters(newscip) );
      SCIP_CALL( SCIPcopyParamSettings(scip, newscip) );
      SCIP_CALL( SCIPsetIntParam(newscip,"presolving/maxrestarts", 0) );
      SCIP_CALL( SCIPsetRealParam(newscip, "limits/time", newtimelimit) );
      SCIP_CALL( SCIPsetLongintParam(newscip, "limits/nodes", -1) );

      newlb = (int) MAX(SCIPceil(scip, dualbound), lb);
      newub = (int) MIN(SCIPfloor(scip, primalbound), ub);

      SCIP_CALL( SCIPcreateModel(newscip, X, Y, &newub, newlb, absmaxX, TRUE, NULL, &success) );
      assert( success );

      /* copy solution */
      newprobdata = SCIPgetProbData(newscip);
      assert( newprobdata != NULL );

      oldprobdata = SCIPgetProbData(scip);
      assert( oldprobdata != NULL );

      oldsol = SCIPgetBestSol(scip);
      SCIP_CALL( SCIPcreateOrigSol(newscip, &newsol, NULL) );

      violatedvars = SCIPprobdataConflictGetViolatedvars(newprobdata);
      isusedvars = SCIPprobdataConflictGetIsusedvars(newprobdata);

      /* copy solution values */
      if ( SCIPisGE(scip, newub, SCIPround(scip, primalbound)) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &usedinequalities, (int) SCIPround(scip, primalbound)) );
         nusedinequalities = 0;

         for (i = 0; i < SCIPprobdataGetNVars(oldprobdata); ++i)
         {
            solvar = SCIPprobdataGetVars(oldprobdata)[i];
            assert( solvar != NULL );

            /* skip unused variables */
            if ( SCIPgetSolVal(scip, oldsol, solvar) < 0.5 )
               continue;

            usedinequalities[nusedinequalities] = i;
            vardata = SCIPvarGetData(solvar);
            assert( vardata != NULL );

            /* set that inequality is used */
            SCIP_CALL( SCIPsetSolVal(newscip, newsol, isusedvars[nusedinequalities], 1.0) );

            /* set cut pattern */
            pattern = SCIPvardataGetConsids(vardata);
            len = SCIPvardataGetNConsids(vardata);
            for (j = 0; j < len; ++j)
            {
               SCIP_CALL( SCIPsetSolVal(newscip, newsol, violatedvars[pattern[j]][nusedinequalities], 1.0) );
            }
            ++nusedinequalities;
         }
         assert( nusedinequalities == SCIPround(scip, primalbound) );

         SCIPfreeBufferArray(scip, &usedinequalities);

         SCIP_CALL( SCIPaddSol(newscip, newsol, &stored) );
         SCIP_CALL( SCIPfreeSol(newscip, &newsol) );
      }

      SCIP_CALL( SCIPsetBoolParam(newscip, "propagating/convexity/enabled", TRUE) );
      SCIP_CALL( SCIPsetBoolParam(newscip, "propagating/intersection/enabled", TRUE) );
      SCIP_CALL( SCIPsolve(newscip) );
      SCIP_CALL( SCIPprintStatistics(newscip, NULL) );

      SCIP_CALL( writeSolutionConflict(newscip, filename, SCIPround(newscip, SCIPgetPrimalbound(newscip)), X, Y, absmaxX) );

      /* free new structures */
      SCIP_CALL( SCIPfreeTransform(newscip) );
      SCIP_CALL( SCIPfree(&newscip) );
   }
   else
   {
      if ( method != METHOD_CG )
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "propagating/convexity/enabled", TRUE) );
         SCIP_CALL( SCIPsetBoolParam(scip, "propagating/intersection/enabled", TRUE) );
      }
      else
      {
         SCIP_CALL( SCIPsetBoolParam(scip, "propagating/convexity/enabled", FALSE) );
         SCIP_CALL( SCIPsetBoolParam(scip, "propagating/intersection/enabled", FALSE) );
         SCIP_CALL( SCIPsetIntParam(scip,"presolving/maxrestarts", 0) );
         SCIP_CALL( SCIPsetSeparating(scip, SCIP_PARAMSETTING_OFF, TRUE) );
      }

      /* use static orbitopes to avoid conflicts with prop_convexity */
      /* SCIP_CALL( SCIPsetBoolParam(scip, "constraints/orbitope/usedynamicprop", FALSE) ); */
      SCIP_CALL( SCIPsolve(scip) );

      /* output statistics */
      SCIPinfoMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );

      if ( method == METHOD_CG )
      {
         SCIP_CALL( writeSolutionCG(scip, filename) );
      }
      else if ( method == METHOD_COMPACT_MIP )
      {
         SCIP_CALL( writeSolutionCompact(scip, filename, ub, X->dimension) );
      }
      else
      {
         SCIP_CALL( writeSolutionConflict(scip, filename, ub, X, Y, absmaxX) );
      }

   }

   return SCIP_OKAY;
}

/** fills existing solution of the problem based on conv(X) based on CDD representation */
SCIP_RETCODE SCIPgetSolutionCompactModelCDD(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   dd_MatrixPtr          facetsconvexhull    /**< facet description of conv(X) */
   )
{
   Datapoints* Y;
   SCIP_VAR*** lhsvars;
   SCIP_VAR** rhsvars;
   SCIP_VAR*** violatedvars;
   SCIP_VAR** isusedvars;
   SCIP_CONS*** linkviolconss;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real eps;
   SCIP_Real maxval;
   int nY;
   int ub;
   int dim;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( facetsconvexhull != NULL );

   lhsvars = SCIPprobdataGetLhsvars(probdata);
   rhsvars = SCIPprobdataGetRhsvars(probdata);
   violatedvars = SCIPprobdataGetViolatedvars(probdata);
   isusedvars = SCIPprobdataGetUsedvars(probdata);
   linkviolconss = SCIPprobdataGetLinkviolconss(probdata);
   Y = SCIPprobdataCompactGetY(probdata);
   nY = Y->ndatapoints;
   ub = SCIPprobdataGetUb(probdata);
   dim = Y->dimension;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* set inequality variables */
   for (i = 0; i < ub; ++i)
   {
      /* find absolute maximum coefficient and possibly rescale inequality */
      maxval = 1.0;
      for (j = 0; j <= dim; ++j)
      {
         if ( SCIPisGT(scip, ABS(getReal(facetsconvexhull->matrix[i][j])), maxval) )
            maxval = ABS(getReal(facetsconvexhull->matrix[i][j]));
      }

      SCIP_CALL( SCIPsetSolVal(scip, sol, rhsvars[i], getReal(facetsconvexhull->matrix[i][0]) / maxval) );

      for (j = 0; j < dim; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, lhsvars[i][j], - getReal(facetsconvexhull->matrix[i][j+1]) / maxval) );
      }
   }

   /* set isusedvars */
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, isusedvars[i], 1.0) );
   }

   /* set violated vars and slack variable of indicator constraints */
   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < ub; ++j)
      {
         val = getReal(facetsconvexhull->matrix[j][0]);

         for (k = 0; k < dim; ++k)
            val += getReal(facetsconvexhull->matrix[j][k+1]) * Y->points[i][k];

         if ( SCIPisLE(scip, val, -eps) )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, violatedvars[i][j], 1.0) );
         }
         else
         {
            var = SCIPgetSlackVarIndicator(linkviolconss[i][j]);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, MAX(val, 0.0) + eps) );
         }
      }
   }

   return SCIP_OKAY;
}

/** fills existing solution of the problem based on explicit solution */
SCIP_RETCODE SCIPgetSolutionCompactModelExplicit(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality */
   int                   ninequalities       /**< number of inequalities encoded in previous data structures */
   )
{
   Datapoints* Y;
   SCIP_VAR*** lhsvars;
   SCIP_VAR** rhsvars;
   SCIP_VAR*** violatedvars;
   SCIP_VAR** isusedvars;
   SCIP_CONS*** linkviolconss;
   SCIP_VAR* var;
   SCIP_Real val;
   SCIP_Real eps;
   int nY;
#ifndef NDEBUG
   int ub;
#endif
   int dim;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( inequalities != NULL );
   assert( separatedpoints != NULL );
   assert( nseparatedpoints != NULL );
   assert( ninequalities > 0 );

   lhsvars = SCIPprobdataGetLhsvars(probdata);
   rhsvars = SCIPprobdataGetRhsvars(probdata);
   violatedvars = SCIPprobdataGetViolatedvars(probdata);
   isusedvars = SCIPprobdataGetUsedvars(probdata);
   linkviolconss = SCIPprobdataGetLinkviolconss(probdata);
   Y = SCIPprobdataCompactGetY(probdata);
   nY = Y->ndatapoints;
#ifndef NDEBUG
   ub = SCIPprobdataGetUb(probdata);
#endif
   dim = Y->dimension;

   assert( ninequalities <= ub );

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* set inequality variables */
   for (i = 0; i < ninequalities; ++i)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, rhsvars[i], inequalities[i][0]) );
      for (j = 0; j < dim; ++j)
      {
         SCIP_CALL( SCIPsetSolVal(scip, sol, lhsvars[i][j], - inequalities[i][j+1]) );
      }
   }

   /* set isusedvars */
   for (i = 0; i < ninequalities; ++i)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, isusedvars[i], 1.0) );
   }

   /* set violated vars and slack variable of indicator constraints */
   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < ninequalities; ++j)
      {
         val = inequalities[j][0];

         for (k = 0; k < dim; ++k)
            val += inequalities[j][k+1] * Y->points[i][k];

         if ( SCIPisLE(scip, val, -eps) )
         {
#ifndef NDEBUG
            for (k = 0; k < nseparatedpoints[j]; ++k)
            {
               if ( separatedpoints[j][k] == i )
                  break;
            }
            assert( k < nseparatedpoints[j] );
#endif
            SCIP_CALL( SCIPsetSolVal(scip, sol, violatedvars[i][j], 1.0) );
         }
         else
         {
#ifndef NDEBUG
            for (k = 0; k < nseparatedpoints[j]; ++k)
            {
               if ( separatedpoints[j][k] == i )
                  break;
            }
            assert( k == nseparatedpoints[j] );
#endif
            var = SCIPgetSlackVarIndicator(linkviolconss[i][j]);
            SCIP_CALL( SCIPsetSolVal(scip, sol, var, MAX(val, 0.0) + eps) );
         }
      }
   }

   return SCIP_OKAY;
}

/** fills existing solution of the problem based on conv(X) based on CDD representation */
SCIP_RETCODE SCIPgetSolutionConflictModelCDD(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   dd_MatrixPtr          facetsconvexhull    /**< facet description of conv(X) */
   )
{
   Datapoints* Y;
   SCIP_VAR*** violatedvars;
   SCIP_VAR** isusedvars;
   SCIP_Real val;
   SCIP_Real eps;
   int nY;
   int ub;
   int dim;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( facetsconvexhull != NULL );

   violatedvars = SCIPprobdataConflictGetViolatedvars(probdata);
   isusedvars = SCIPprobdataConflictGetIsusedvars(probdata);
   Y = SCIPprobdataConflictGetY(probdata);
   nY = Y->ndatapoints;
   ub = SCIPprobdataConflictGetUb(probdata);
   dim = Y->dimension;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* set isusedvars */
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, isusedvars[i], 1.0) );
   }

   /* set violated vars and slack variable of indicator constraints */
   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < ub; ++j)
      {
         val = getReal(facetsconvexhull->matrix[j][0]);

         for (k = 0; k < dim; ++k)
            val += getReal(facetsconvexhull->matrix[j][k+1]) * Y->points[i][k];

         if ( SCIPisLE(scip, val, -eps) )
         {
            SCIP_CALL( SCIPsetSolVal(scip, sol, violatedvars[i][j], 1.0) );
         }
      }
   }

   return SCIP_OKAY;
}

/** fills existing solution of the problem based on explicit solution */
SCIP_RETCODE SCIPgetSolutionConflictModelExplicit(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_SOL*             sol,                /**< solution to be filled */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality */
   int                   ninequalities       /**< number of inequalities encoded in previous data structures */
   )
{
   Datapoints* Y;
   SCIP_VAR*** violatedvars;
   SCIP_VAR** isusedvars;
   SCIP_Real val;
   SCIP_Real eps;
   int nY;
#ifndef NDEBUG
   int ub;
#endif
   int dim;
   int i;
   int j;
   int k;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( inequalities != NULL );
   assert( separatedpoints != NULL );
   assert( nseparatedpoints != NULL );
   assert( ninequalities > 0 );

   violatedvars = SCIPprobdataConflictGetViolatedvars(probdata);
   isusedvars = SCIPprobdataConflictGetIsusedvars(probdata);
   Y = SCIPprobdataConflictGetY(probdata);
   nY = Y->ndatapoints;
#ifndef NDEBUG
   ub = SCIPprobdataConflictGetUb(probdata);
#endif
   dim = Y->dimension;

   assert( ninequalities <= ub );

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* set isusedvars */
   for (i = 0; i < ninequalities; ++i)
   {
      SCIP_CALL( SCIPsetSolVal(scip, sol, isusedvars[i], 1.0) );
   }

   /* set violated vars and slack variable of indicator constraints */
   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < ninequalities; ++j)
      {
         val = inequalities[j][0];

         for (k = 0; k < dim; ++k)
            val += inequalities[j][k+1] * Y->points[i][k];

         if ( SCIPisLE(scip, val, -eps) )
         {
#ifndef NDEBUG
            for (k = 0; k < nseparatedpoints[j]; ++k)
            {
               if ( separatedpoints[j][k] == i )
                  break;
            }
            assert( k < nseparatedpoints[j] );
#endif
            SCIP_CALL( SCIPsetSolVal(scip, sol, violatedvars[i][j], 1.0) );
         }
      }
   }

   return SCIP_OKAY;
}
