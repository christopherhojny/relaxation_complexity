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
 * @brief  functions for computing maximally separable sets
 * @author Christopher Hojny
 */

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"

#include "cons_samediff.h"
#include "maximal_separation.h"
#include "typedefs.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"


/** creates variables of separation problem */
static
SCIP_RETCODE createVars(
   SCIP*                 subscip,            /**< SCIP pointer of separation model */
   SCIP_VAR**            ineqvars,           /**< variables modeling separating inequality */
   SCIP_VAR**            isviolvars,         /**< variables modeling violation */
   int                   dimension,          /**< dimension of data points */
   int                   nY,                 /**< number of points to be separated */
   SCIP_Real*            objvals,            /**< objective values of points to be separated (or NULL for all ones) */
   int                   absmaxX             /**< maximum absolute value of a coordinate in feasible points */
   )
{
   char name[SCIP_MAXSTRLEN];
   int bound;
   int i;

   assert( subscip != NULL );
   assert( ineqvars != NULL );
   assert( isviolvars != NULL );
   assert( dimension > 0 );
   assert( nY > 0 );

   /* create variables encoding separating inequality */
   for (i = 0; i < dimension; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "a_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &ineqvars[i], name, -1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, ineqvars[i]) );
   }
   bound = dimension * absmaxX;
   SCIP_CALL( SCIPcreateVar(subscip, &ineqvars[dimension], "b", -bound, bound, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, ineqvars[dimension]) );

   /* create variables indicating whether point is separated */
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "viol_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &isviolvars[i], name, 0.0, 1.0, objvals == NULL ? 1.0 : objvals[i],
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, isviolvars[i]) );
   }

   return SCIP_OKAY;
}

/** releases variables of separation model */
static
SCIP_RETCODE releaseVars(
   SCIP*                 subscip,            /**< SCIP pointer of separation model */
   SCIP_VAR**            ineqvars,           /**< variables modeling separating inequality */
   SCIP_VAR**            isviolvars,         /**< variables modeling violation */
   int                   dimension,          /**< dimension of data points */
   int                   nY                  /**< number of points to be separated */
   )
{
   int i;

   assert( subscip != NULL );
   assert( ineqvars != NULL );
   assert( isviolvars != NULL );
   assert( dimension > 0 );
   assert( nY > 0 );

   for (i = 0; i < nY ; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &isviolvars[i]) );
   }

   for (i = 0; i <= dimension ; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &ineqvars[i]) );
   }

   return SCIP_OKAY;
}

/** creates constraints of separation model */
static
SCIP_RETCODE createConss(
   SCIP*                 subscip,            /**< SCIP pointer of separation model */
   SCIP_CONS**           validconss,         /**< array for constraints enforcing validity of separating inequality */
   SCIP_CONS**           sepaconss,          /**< array for constraints enforcing separation property */
   SCIP_VAR**            ineqvars,           /**< variables modeling separating inequality */
   SCIP_VAR**            isviolvars,         /**< variables modeling violation */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   nX,                 /**< number of feasible points */
   int                   nY,                 /**< number of points to be separated */
   int                   dimension,          /**< dimension of data points */
   SCIP_Real             eps                 /**< epsilon for classifying violation */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real* tmpvals;
   int i;
   int j;

   assert( subscip != NULL );
   assert( validconss != NULL );
   assert( sepaconss != NULL );
   assert( ineqvars != NULL );
   assert( isviolvars != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( dimension > 0);

   SCIP_CALL( SCIPallocBufferArray(subscip, &tmpvals, dimension + 1) );
   tmpvals[dimension] = -1.0;

   /* create constraints enforcing validity */
   for (i = 0; i < nX; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "validcons_%d", i);

      for (j = 0; j < dimension; ++j)
         tmpvals[j] = X->points[i][j];

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &validconss[i], name,
            dimension + 1, ineqvars, tmpvals, -SCIPinfinity(subscip), 0.0) );
      SCIP_CALL( SCIPaddCons(subscip, validconss[i]) );
   }

   /* create constraints modeling violation */
   tmpvals[dimension] = 1.0;
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sepacons_%d", i);

      for (j = 0; j < dimension; ++j)
         tmpvals[j] = -Y->points[i][j];

      SCIP_CALL( SCIPcreateConsBasicIndicator(subscip, &sepaconss[i], name,
            isviolvars[i], dimension + 1, ineqvars, tmpvals, -eps) );
      SCIP_CALL( SCIPaddCons(subscip, sepaconss[i]) );
   }

   SCIPfreeBufferArray(subscip, &tmpvals);

   return SCIP_OKAY;
}

/** releases constraints of separation model */
static
SCIP_RETCODE releaseConss(
   SCIP*                 subscip,            /**< SCIP pointer of separation model */
   SCIP_CONS**           validconss,         /**< array for constraints enforcing validity of separating inequality */
   SCIP_CONS**           sepaconss,          /**< array for constraints enforcing separation property */
   int                   nX,                 /**< number of feasible points */
   int                   nY,                 /**< number of points to be separated */
   int                   dimension           /**< dimension of data points */
   )
{
   int i;

   assert( subscip != NULL );
   assert( validconss != NULL );
   assert( sepaconss != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( dimension > 0);

   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &sepaconss[i]) );
   }
   for (i = 0; i < nX; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &validconss[i]) );
   }

   return SCIP_OKAY;
}

/** applies fixings that require a point to (not) be separated */
static
SCIP_RETCODE applyFixings(
   SCIP*                 subscip,            /**< SCIP pointer of separation model */
   SCIP_VAR**            isviolvars,         /**< variables modeling violation */
   int*                  fixedY,             /**< array of points in Y whose separation status is fixed
                                                  (or NULL if not needed)*/
   SCIP_Real*            statusfixedY,       /**< separation status of points in fixedY (or NULL if not needed):
                                                0 - not separated, 1 - separated */
   int                   nfixedY             /**< number of points in fixedY */
   )
{
   int i;

   assert( subscip != NULL );
   assert( isviolvars != NULL );
   assert( fixedY != NULL );
   assert( statusfixedY != NULL );
   assert( nfixedY > 0 );

   for (i = 0; i < nfixedY; ++i)
   {
      assert( SCIPisEQ(subscip, statusfixedY[i], 1.0) || SCIPisEQ(subscip, statusfixedY[i], 0.0) );

      if ( SCIPisEQ(subscip, statusfixedY[i], 1.0) )
      {
         SCIP_CALL( SCIPchgVarLbGlobal(subscip, isviolvars[fixedY[i]], 1.0) );
      }
      else
      {
         SCIP_CALL( SCIPchgVarUbGlobal(subscip, isviolvars[fixedY[i]], 0.0) );
      }
   }

   return SCIP_OKAY;
}

/** add branching decisions constraints to the sub SCIP */
static
SCIP_RETCODE addBranchingDecisionConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array for violation vars */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for branching data (or NULL) */
   )
{
   SCIP_CONS** conss;
   SCIP_CONS* cons;
   int nconss;
   int id1;
   int id2;
   CONSTYPE type;
   SCIP_Real vbdcoef;
   SCIP_Real lhs;
   SCIP_Real rhs;
   int c;

   if ( conshdlr == NULL )
      return SCIP_OKAY;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( conshdlr != NULL );

   /* get branching decisions */
   conss = SCIPconshdlrGetConss(conshdlr);
   nconss = SCIPconshdlrGetNConss(conshdlr);

   /* loop over all branching decision constraints and apply the branching decision if the corresponding constraint is
    * active
    */
   for (c = 0; c < nconss; ++c)
   {
      cons = conss[c];

      /* ignore inactive constraints (not laying on the path from the current node to the root) */
      if ( !SCIPconsIsActive(cons) )
         continue;

      /* get the IDs of patterns that have been used for branching and the branching type (SAME or DIFFER) */
      id1 = SCIPgetItemid1Samediff(scip, cons);
      id2 = SCIPgetItemid2Samediff(scip, cons);
      type = SCIPgetTypeSamediff(scip, cons);

      /* depending on the branching type select the correct left and right hand side for the linear constraint which
       * enforces this branching decision in the pricing problem MIP
       */
      if ( type == SAME )
      {
         lhs = 0.0;
         rhs = 0.0;
         vbdcoef = -1.0;
      }
      else if ( type == DIFFER )
      {
         lhs = -SCIPinfinity(scip);
         rhs = 1.0;
         vbdcoef = 1.0;
      }
      else
      {
         SCIPerrorMessage("unknow constraint type <%d>\n", type);
         return SCIP_INVALIDDATA;
      }

      /* add linear (in that case a variable bound) constraint to pricing MIP depending on the branching type:
       *
       * - branching type SAME:  x1 = x2 <=> x1 - x2 = 0 <=> 0 <= x1 - x2 <= 0
       *
       * - branching type DIFFER:  x1 + x2 <= 1 <=> -inf <= x1 + x2 <= 1
       *
       * note a setppc constraint would be sufficient and even better suitable for such kind of constraint
       */
      SCIP_CALL( SCIPcreateConsBasicVarbound(subscip, &cons, SCIPconsGetName(conss[c]),
            vars[id1], vars[id2], vbdcoef, lhs, rhs) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   return SCIP_OKAY;
}

/** finds an inequality separating as many infeasible points as possible */
SCIP_RETCODE SCIPfindMaximalSeparatingInequality(
   SCIP*                 scip,               /**< main SCIP pointer */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   nX,                 /**< number of data points encoded in X */
   int                   nY,                 /**< number of data points encoded in Y */
   SCIP_Real*            objvals,            /**< objective values of points to be separated (or NULL for all ones) */
   int                   dimension,          /**< dimension of data points */
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   int*                  fixedY,             /**< array of points in Y whose separation status is fixed
                                                  (or NULL if not needed)*/
   SCIP_Real*            statusfixedY,       /**< separation status of points in fixedY (or NULL if not needed):
                                                0 - not separated, 1 - separated */
   int                   nfixedY,            /**< number of points in fixedY */
   int*                  separatedY,         /**< array of points separated in solution (needs to be allocated) */
   int*                  nseparatedY,        /**< pointer to store number of separated points */
   SCIP_Real*            sepainequality,     /**< array to store separated inequality (needs to be allocated) */
   SCIP_CONSHDLR*        conshdlr,           /**< conshdlr for branching decisions (or NULL if not needed) */
   SCIP_Longint          nodelimit,          /**< node limit for solving the separation problem (or -1 if no limit) */
   SCIP_Real             timelimit,          /**< time limit for finding maximumally violated inequality (of -1 if no limit) */
   SCIP_Bool*            success             /**< pointer to store whether a separating inequality could be found */
   )
{
   SCIP* subscip;
   SCIP_SOL* sol;
   SCIP_VAR** ineqvars;
   SCIP_VAR** isviolvars;
   SCIP_CONS** validconss;
   SCIP_CONS** sepaconss;
   SCIP_Real eps;
   int i;
   int j;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( dimension > 0 );
   assert( fixedY != NULL || nfixedY == 0 );
   assert( statusfixedY != NULL || nfixedY == 0 );
   assert( separatedY != NULL );
   assert( nseparatedY != NULL );
   assert( sepainequality != NULL );
   assert( success != NULL );

   *success = TRUE;
   *nseparatedY = 0;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* basic setup */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcreateProbBasic(subscip, "separatingInequality") );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* maximize number of separated points */
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* load plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &ineqvars, dimension + 1) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &isviolvars, nY) );
   SCIP_CALL( createVars(subscip, ineqvars, isviolvars, dimension, nY, objvals, absmaxX) );

   /* create constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &validconss, nX) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &sepaconss, nY) );
   SCIP_CALL( createConss(subscip, validconss, sepaconss, ineqvars, isviolvars,
         X, Y, nX, nY, dimension, eps) );

   /* apply fixings */
   if ( nfixedY > 0 )
   {
      SCIP_CALL( applyFixings(subscip, isviolvars, fixedY, statusfixedY, nfixedY) );
   }

   /* incorporate decisions encoded via samediff conss */
   SCIP_CALL( addBranchingDecisionConss(scip, subscip, isviolvars, conshdlr) );

   /* solve problem */
   if ( nodelimit >= 0 )
   {
      SCIP_CALL( SCIPsetLongintParam(subscip, "limits/nodes", nodelimit) );
   }
   if ( SCIPisGE(scip, timelimit, 0.0) )
   {
      SCIP_CALL( SCIPsetRealParam(subscip, "limits/time", timelimit) );
   }
   SCIP_CALL( SCIPsolve(subscip) );

   /* terminate if the user has interupted the pricing problem */
   if ( SCIPpressedCtrlC(subscip) )
   {
      SCIPdebugMsg(scip, "no separating inequality could be found due to user interrupt \n");
      *success = FALSE;

      goto FREEDATA;
   }

   /* terminate if no violated inequality has been found */
   if ( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD )
   {
      SCIPdebugMsg(scip, "no separating inequality could be found\n");
      *success = FALSE;

      goto FREEDATA;
   }

   /* extract solution */
   sol = SCIPgetBestSol(subscip);
   sepainequality[0] = SCIPgetSolVal(subscip, sol, ineqvars[dimension]);
   for (i = 0; i < dimension; ++i)
      sepainequality[i + 1] = -SCIPgetSolVal(subscip, sol, ineqvars[i]);

   for (i = 0; i < nY; ++i)
   {
      SCIP_Real viol;

      viol = sepainequality[0];
      for (j = 0; j < dimension; ++j)
         viol += Y->points[i][j] * sepainequality[j + 1];

      if ( SCIPisLE(subscip, viol, -eps) )
      {
         separatedY[*nseparatedY] = i;
         *nseparatedY += 1;
      }
   }

   FREEDATA:
   SCIP_CALL( releaseConss(subscip, validconss, sepaconss, nX, nY, dimension) );
   SCIPfreeBlockMemoryArrayNull(subscip, &sepaconss, nY);
   SCIPfreeBlockMemoryArrayNull(subscip, &validconss, nX);

   SCIP_CALL( releaseVars(subscip, ineqvars, isviolvars, dimension, nY) );
   SCIPfreeBlockMemoryArrayNull(subscip, &isviolvars, nY);
   SCIPfreeBlockMemoryArrayNull(subscip, &ineqvars, dimension + 1);

   SCIP_CALL( SCIPfreeTransform(subscip) );
   SCIP_CALL( SCIPfree(&subscip) );


   return SCIP_OKAY;
}

/** finds a greedy solution for RC */
SCIP_RETCODE heurGreedy(
   SCIP*                 scip,               /**< SCIP pointer */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   nX,                 /**< number of data points encoded in X */
   int                   nY,                 /**< number of data points encoded in Y */
   int                   dimension,          /**< dimension of data points */
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   int                   ub,                 /**< upper bound on RC */
   int                   greedysort,         /**< sorting of infeasible points used by heuristic */
   int                   method,             /**< method used to compute RC */
   SCIP_Real**           inequalities,       /**< (ub x dimension+1) array to store inequalities */
   int**                 separatedpoints,    /**< (ub x nY) array to store list of separated points */
   int*                  nseparatedpoints,   /**< ub-dimensional array to store for each inequality
                                                the number of separated points */
   int*                  ninequalities,      /**< pointer to store number of inequalities in solution */
   dd_MatrixPtr          facetsconvexhull,   /**< facets of conv(X) (or NULL if not needed) */
   SCIP_Real             timelimit           /**< time limit for greedy computation */
   )
{
   SCIP_Real** heurinequalities;
   int** heurseparatedpoints;
   SCIP_Real* objvals;
   SCIP_Longint nodelimit;
   int* sorting;
   int* tmpseparated;
   int ntmpseparated;
   int i;
   int j;
   SCIP_Bool success;
   int k = 0;
   int heurnseparatedpoints = 0;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( dimension > 0 );
   assert( absmaxX > 0 );
   assert( 0 <= greedysort && greedysort <= 1 );
   assert( SCIPisGE(scip, timelimit, 0.0) );

   *ninequalities = -1;

   SCIP_CALL( SCIPallocBufferArray(scip, &sorting, nY) );
   for (i = 0; i < nY; ++i)
      sorting[i] = i;

   if ( greedysort == GREEDY_SORT_DISTANCE )
   {
      SCIP_Real* latticedistances;
      int c;

      /* cannot compute greedy solution */
      if ( facetsconvexhull == NULL )
      {
         SCIPinfoMessage(scip, NULL, "   cannot compute greedy solution with distance strategy: convex hull not available\n");
         SCIPfreeBufferArray(scip, &sorting);

         return SCIP_OKAY;
      }

      SCIP_CALL( SCIPallocBufferArray(scip, &latticedistances, nY) );
      for (i = 0; i < nY; ++i)
      {
         latticedistances[i] = SCIP_REAL_MAX;

         for (j = 0; j < facetsconvexhull->rowsize; ++j)
         {
            SCIP_Real val;

            /* compute lattice distance to facet */
            val = - getReal(facetsconvexhull->matrix[j][0]);
            for (c = 1; c < facetsconvexhull->colsize; ++c)
               val -= getReal(facetsconvexhull->matrix[j][c]) * Y->points[i][c - 1];

            /* if we have found a facet with smaller distance */
            if ( SCIPisGT(scip, val, 0.0) && SCIPisLT(scip, val, latticedistances[i]) )
               latticedistances[i] = val;
         }
      }

      SCIPsortRealInt(latticedistances, sorting, nY);

      SCIPfreeBufferArray(scip, &latticedistances);
   }

   /* iteratively compute a maximum size set of points of not yet separated points */
   SCIP_CALL( SCIPallocBufferArray(scip, &heurinequalities, ub) );
   SCIP_CALL( SCIPallocBufferArray(scip, &heurseparatedpoints, ub) );
   SCIP_CALL( SCIPallocBufferArray(scip, &objvals, nY) );
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpseparated, nY) );
   for (i = 0; i < nY; ++i)
      objvals[i] = 1.0;

   SCIP_CALL( SCIPgetLongintParam(scip, "rc/heurnodelimit", &nodelimit) );

   while ( k < ub && heurnseparatedpoints < nY )
   {
      SCIP_CALL( SCIPallocBufferArray(scip, &heurinequalities[k], dimension + 1) );

      SCIP_CALL( SCIPfindMaximalSeparatingInequality(scip, X, Y, nX, nY, objvals, dimension, absmaxX,
            NULL, NULL, 0, tmpseparated, &ntmpseparated, heurinequalities[k], NULL, nodelimit, timelimit, &success) );

      if ( ! success )
         break;
      assert( ntmpseparated > 0 );

      SCIP_CALL( SCIPallocBufferArray(scip, &heurseparatedpoints[k], ntmpseparated + 1) );

      /* identify the points separated by the k-th inequality
       *
       * the first entry of separatedpoints[k] stores how many points are separated by the inequality
       */
      for (i = 0; i < ntmpseparated; ++i)
      {
         if ( objvals[tmpseparated[i]] > 0.5 )
         {
            objvals[tmpseparated[i]] = 0.0;
            ++heurnseparatedpoints;
         }
      }
      heurseparatedpoints[k][0] = ntmpseparated;
      for (i = 0; i < ntmpseparated; ++i)
         heurseparatedpoints[k][i + 1] = tmpseparated[i];
      k++;
   }

   /* possibly store a solution */
   if ( k < ub && heurnseparatedpoints == nY )
   {
      SCIPdebugMsg(scip, "found a solution with %d inequalities.\n", k);

      for (i = 0; i < k; ++i)
      {
         for (j = 0; j <= dimension; ++j)
            inequalities[i][j] = heurinequalities[i][j];
      }
      *ninequalities = k;

      for (i = 0; i < k; ++i)
      {
         for (j = 0; j < heurseparatedpoints[i][0]; ++j)
            separatedpoints[i][j] = heurseparatedpoints[i][j+1];
         nseparatedpoints[i] = heurseparatedpoints[i][0];
      }
   }

   /* free memory */
   if ( ! success )
   {
      SCIPfreeBufferArray(scip, &heurinequalities[k]);
   }
   --k;
   while ( k >= 0 )
   {
      SCIPfreeBufferArray(scip, &heurseparatedpoints[k]);
      SCIPfreeBufferArray(scip, &heurinequalities[k--]);
   }

   SCIPfreeBufferArray(scip, &tmpseparated);
   SCIPfreeBufferArray(scip, &objvals);
   SCIPfreeBufferArray(scip, &heurseparatedpoints);
   SCIPfreeBufferArray(scip, &heurinequalities);
   SCIPfreeBufferArray(scip, &sorting);

   return SCIP_OKAY;
}
