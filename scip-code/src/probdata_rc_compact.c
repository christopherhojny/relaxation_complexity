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

/**@file   probdata_rc_compact.c
 * @brief  Problem data for computing RC using a compact model
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "compute_symmetries.h"
#include "hiding_sets.h"
#include "datapoints.h"
#include "probdata_rc_compact.h"
#include "problem_rc.h"
#include "vardata_compact.h"

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"
#include "math.h"
#include "convex_hull.h"

#include "scip/scip.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_orbitope.h"
#include "scip/cons_setppc.h"
#include "scip/cons_symresack.h"
#include "scip/cons_varbound.h"

/** @brief Problem data which is accessible in all places
 *
 * This problem data is used to store the input of the problem, all variables which are created, and all
 * constraints.
 */
struct SCIP_ProbData
{
   /* data of the problem */
   Datapoints*           X;                  /**< pointer to data points of X */
   Datapoints*           Y;                  /**< pointer to data points of Y */
   int                   nX;                 /**< number of data points encoded in X */
   int                   nY;                 /**< number of data points encoded in Y */
   int                   ub;                 /**< upper bound on RC(X,Y) */
   int                   lb;                 /**< lower bound on RC(X,Y) */
   int                   dimension;          /**< dimension of data points */
   int                   absmaxX;            /**< maximum absolute entry of a coordinate in X */
   int**                 perms;              /**< permutation symmetries of Y w.r.t. X */
   int                   nperms;             /**< number of permutations stored in perms */
   int                   nmaxperms;          /**< maximum number of permutations that fit into perms */

   /* data of the model */
   SCIP_VAR***           lhsvars;            /**< left-hand side variables of relaxation */
   SCIP_VAR**            rhsvars;            /**< right-hand side variables of relaxation */
   SCIP_VAR***           violatedvars;       /**< variables indicating whether constraint is violated */
   SCIP_VAR**            isusedvars;         /**< variables indicating whether constraint is used */
   SCIP_CONS***          validconss;         /**< constraints enforcing relaxation to be valid */
   SCIP_CONS**           violationconss;     /**< constraints ensuring each infeasible point violated a constraint */
   SCIP_CONS***          linkviolconss;      /**< constraints to link violatedvars with constaints */
   SCIP_CONS***          linkusedconss;      /**< constraints to link isusedvars with constraints */
   SCIP_CONS*            objbound;           /**< constraint modeling lower bound on objective */

   /* advanced constraints */
   SCIP_CONS*            orbitopecons;       /**< constraint to handle permutations of infeasible point classes */
   SCIP_CONS**           symresackconss;     /**< constraints to handle permutations of Y */
   SCIP_CONS**           sortusedconss;      /**< constraints to enforce sorting on isusedvars */
   SCIP_CONS**           sortdifflhsconss;   /**< constraints to enforce sorting between different constraints */
   SCIP_CONS***          uppervbdconss;      /**< constraints to force upper bound on lhsvars to 0 if inequality not used */
   SCIP_CONS***          lowervbdconss;      /**< constraints to force lower bound on lhsvars to 0 if inequality not used */
   SCIP_CONS**           rhsvbdconss;        /**< constraints to force rhs to maximum value if inequality is not used */
   SCIP_CONS**           hidingsetcuts;      /**< hiding set based cuts */
   int                   nhidingsetcuts;     /**< number of hiding set cuts stored in hidingsetcuts */
};


/**@name Local methods
 *
 * @{
 */


/** creates problem data */
static
SCIP_RETCODE probdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata,           /**< pointer to problem data */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   ub,                 /**< upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   int**                 perms,              /**< permutation symmetries of Y w.r.t. X */
   int                   nperms,             /**< number of permutations stored in perms */
   int                   nmaxperms,          /**< maximum number of permutations that fit into perms */
   SCIP_VAR***           lhsvars,            /**< left-hand side variables of relaxation (ub x dimension) */
   SCIP_VAR**            rhsvars,            /**< right-hand side variables of relaxation */
   SCIP_VAR***           violatedvars,       /**< variables indicating whether constraint is violated (Y x ub)*/
   SCIP_VAR**            isusedvars,         /**< variables indicating whether constraint is used */
   SCIP_CONS***          validconss,         /**< constraints enforcing relaxation to be valid */
   SCIP_CONS**           violationconss,     /**< constraints ensuring each infeasible point violated a constraint */
   SCIP_CONS***          linkviolconss,      /**< constraints to link violatedvars with constaints (Y x ub) */
   SCIP_CONS***          linkusedconss,      /**< constraints to link isusedvars with constraints (Y x ub) */
   SCIP_CONS*            objbound,           /**< constraint modeling lower bound on objective */
   SCIP_CONS*            orbitopecons,       /**< constraint to handle permutations of infeasible point classes */
   SCIP_CONS**           symresackconss,     /**< constraints to handle permutations of Y */
   SCIP_CONS**           sortusedconss,      /**< constraints to enforce sorting on isusedvars */
   SCIP_CONS**           sortdifflhsconss,   /**< constraints to enforce sorting between different constraints */
   SCIP_CONS***          uppervbdconss,      /**< constraints to force upper bound on lhsvars to 0 if inequality not used */
   SCIP_CONS***          lowervbdconss,      /**< constraints to force lower bound on lhsvars to 0 if inequality not used */
   SCIP_CONS**           rhsvbdconss,        /**< constraints to force rhs to maximum value if inequality is not used */
   SCIP_CONS**           hidingsetcuts,      /**< hiding set based cuts */
   int                   nhidingsetcuts      /**< number of hiding set cuts stored in hidingsetcuts */
   )
{
   int nX;
   int nY;
   int dimension;
   int i;
   SCIP_Bool decision;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub > 0 );

   /* allocate memory */
   SCIP_CALL( SCIPallocBlockMemory(scip, probdata) );

   nX = X->ndatapoints;
   nY = Y->ndatapoints;
   dimension = X->dimension;

   (*probdata)->X = X;
   (*probdata)->Y = Y;
   (*probdata)->nX = nX;
   (*probdata)->nY = nY;
   (*probdata)->dimension = dimension;
   (*probdata)->ub = ub;
   (*probdata)->lb = lb;
   (*probdata)->absmaxX = absmaxX;
   (*probdata)->nhidingsetcuts = nhidingsetcuts;
   (*probdata)->nperms = nperms;
   (*probdata)->nmaxperms = nmaxperms;

   (*probdata)->objbound = objbound;

   /* possible copy variable arrays */
   if ( lhsvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->lhsvars, lhsvars, ub) );
      for (i = 0; i < ub; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->lhsvars[i], lhsvars[i], dimension) );
      }
   }
   else
      (*probdata)->lhsvars = NULL;

   if ( rhsvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rhsvars, rhsvars, ub) );
   }
   else
      (*probdata)->rhsvars = NULL;

   if ( violatedvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->violatedvars, violatedvars, nY) );
      for (i = 0; i < nY; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->violatedvars[i], violatedvars[i], ub) );
      }
   }
   else
      (*probdata)->violatedvars = NULL;

   if ( isusedvars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->isusedvars, isusedvars, ub) );
   }
   else
      (*probdata)->isusedvars = NULL;

   /* duplicate constraint arrays */
   if ( validconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->validconss, validconss, ub) );
      for (i = 0; i < ub; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->validconss[i], validconss[i], nX) );
      }
   }
   else
      (*probdata)->validconss = NULL;

   if ( violationconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->violationconss, violationconss, nY) );
   }
   else
      (*probdata)->violationconss = NULL;

   if ( linkviolconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->linkviolconss, linkviolconss, nY) );
      for (i = 0; i < nY; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->linkviolconss[i], linkviolconss[i], ub) );
      }
   }
   else
      linkviolconss = NULL;

   if ( linkusedconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->linkusedconss, linkusedconss, nY) );
      for (i = 0; i < nY; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->linkusedconss[i], linkusedconss[i], ub) );
      }
   }
   else
      linkusedconss = NULL;

   if ( uppervbdconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->uppervbdconss, uppervbdconss, ub) );
      for (i = 0; i < ub; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->uppervbdconss[i], uppervbdconss[i], dimension) );
      }
   }
   else
      uppervbdconss = NULL;
   if ( lowervbdconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->lowervbdconss, lowervbdconss, ub) );
      for (i = 0; i < ub; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->lowervbdconss[i], lowervbdconss[i], dimension) );
      }
   }
   else
      lowervbdconss = NULL;
   if ( rhsvbdconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->rhsvbdconss, rhsvbdconss, ub) );
   }
   else
      rhsvbdconss = NULL;

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      if ( sortusedconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->sortusedconss, sortusedconss, ub - 1) );
      }
      else
         (*probdata)->sortusedconss = NULL;

      if ( sortdifflhsconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->sortdifflhsconss, sortdifflhsconss, ub - 1) );
      }
      else
         (*probdata)->sortdifflhsconss = NULL;

      if ( symresackconss != NULL )
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->symresackconss, symresackconss, nperms) );
      }
      else
         (*probdata)->symresackconss = NULL;
   }
   (*probdata)->orbitopecons = orbitopecons;

   if ( hidingsetcuts != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->hidingsetcuts, hidingsetcuts, nhidingsetcuts) );
   }
   else
      (*probdata)->hidingsetcuts = NULL;

   if ( perms != NULL )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*probdata)->perms, (*probdata)->nmaxperms) );
      for (i = 0; i < (*probdata)->nperms; ++i)
      {
         SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->perms[i], perms[i], (*probdata)->nY) );
      }
   }
   else
      (*probdata)->perms = NULL;

   return SCIP_OKAY;
}

/** frees the memory of the given problem data */
static
SCIP_RETCODE probdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA**       probdata            /**< pointer to problem data */
   )
{
   int i;
   int j;
   int nX;
   int nY;
   int ub;
   int dimension;
   SCIP_Bool decision;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->X != NULL );
   assert( (*probdata)->Y != NULL );
   assert( (*probdata)->ub > 0 );

   nX = (*probdata)->nX;
   nY = (*probdata)->nY;
   ub = (*probdata)->ub;
   dimension = (*probdata)->dimension;

   /* release all variables */
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->isusedvars[i]) );
   }
   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < ub; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->violatedvars[i][j]) );
      }
   }
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->rhsvars[i]) );
   }
   for (i = 0; i < ub; ++i)
   {
      for (j = 0; j < (*probdata)->dimension; ++j)
      {
         SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->lhsvars[i][j]) );
      }
   }

   /* release all constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->objbound) );

   for (i = nY - 1; i >= 0;--i)
   {
      for (j = 0; j < ub; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linkusedconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->linkusedconss[i], ub);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->linkusedconss, nY);

   for (i = nY - 1; i >= 0;--i)
   {
      for (j = 0; j < ub; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->linkviolconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->linkviolconss[i], ub);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->linkviolconss, nY);

   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->violationconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->violationconss, nY);

   for (i = ub - 1; i >= 0; --i)
   {
      for (j = 0; j < nX; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->validconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->validconss[i], nX);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->validconss, ub);

   for (i = ub - 1; i >= 0; --i)
   {
      for (j = 0; j < dimension; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->uppervbdconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->uppervbdconss[i], dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->uppervbdconss, ub);
   for (i = ub - 1; i >= 0; --i)
   {
      for (j = 0; j < dimension; ++j)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->lowervbdconss[i][j]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->lowervbdconss[i], dimension);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->lowervbdconss, ub);
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->rhsvbdconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rhsvbdconss, ub);


   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "rc/useorbitope", &decision) );

      if ( decision )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->orbitopecons) );

         SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision) );
         if ( decision )
         {
            for (i = 0; i < (*probdata)->nperms; ++i)
            {
               SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->symresackconss[i]) );
            }
            SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->symresackconss, (*probdata)->nperms);
         }
      }
      else
      {
         for (i = 0; i < ub - 1; ++i)
         {
            SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->sortdifflhsconss[i]) );
         }
         SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->sortdifflhsconss, ub - 1);
      }
      for (i = 0; i < ub - 1; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->sortusedconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->sortusedconss, ub - 1);
   }

   for (i = 0; i < (*probdata)->nhidingsetcuts; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->hidingsetcuts[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->hidingsetcuts, (*probdata)->nhidingsetcuts);

   /* free memory of arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->isusedvars, ub);
   for (i = nY - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->violatedvars[i], ub);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->violatedvars, nY);
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->rhsvars, ub);
   for (i = ub - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->lhsvars[i], (*probdata)->dimension);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->lhsvars, ub);
   for (i = (*probdata)->nperms - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->perms[i], (*probdata)->nY);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->perms, (*probdata)->nmaxperms);

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}

/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   int                   absmaxX             /**< maximum absolute value of a coordinate in X */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real varbd;
   int nY;
   int dimension;
   int ub;
   int i;
   int j;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->nX > 0 );
   assert( probdata->nY > 0 );
   assert( probdata->dimension > 0 );
   assert( probdata->ub > 0 );

   nY = probdata->nY;
   dimension = probdata->dimension;
   ub = probdata->ub;

   varbd = (SCIP_Real) (dimension * absmaxX);

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->lhsvars, ub) );
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->lhsvars[i], dimension) );

      for (j = 0; j < dimension; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "a_%d_%d", i, j);

         SCIP_CALL( SCIPcreateVar(scip, &probdata->lhsvars[i][j], name, -1.0, 1.0, 0.0,
               SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(scip, probdata->lhsvars[i][j]) );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->rhsvars, ub) );
   for (i = 0; i < ub; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->rhsvars[i], name, - varbd, varbd, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->rhsvars[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->violatedvars, nY) );
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->violatedvars[i], ub) );

      for (j = 0; j < ub; ++j)
      {
         SCIP_VARDATA* vardata;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "isviol_%d_%d", i, j);
         SCIP_CALL( SCIPvardataCreateCompact(scip, &vardata, j) );

         SCIP_CALL( SCIPcreateVarIsviolatedCompact(scip, &probdata->violatedvars[i][j], name, 0.0, vardata) );
         SCIP_CALL( SCIPaddVar(scip, probdata->violatedvars[i][j]) );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->isusedvars, ub) );
   for (i = 0; i < ub; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "isused_%d", i);

      SCIP_CALL( SCIPcreateVar(scip, &probdata->isusedvars[i], name, 0.0, 1.0, 1.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(scip, probdata->isusedvars[i]) );
   }

   return SCIP_OKAY;
}

/** creates the constraints of the problem */
static
SCIP_RETCODE SCIPcreateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nX;
   int nY;
   int dimension;
   int ub;
   int lb;
   int i;
   int j;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   SCIP_Bool decision;
   int maxnvars;
   SCIP_Real eps;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->X != NULL );
   assert( probdata->Y != NULL );
   assert( probdata->X->points != NULL );
   assert( probdata->Y->points != NULL );
   assert( probdata->nX > 0 );
   assert( probdata->nY > 0 );
   assert( probdata->dimension > 0 );
   assert( probdata->ub > 0 );

   nX = probdata->nX;
   nY = probdata->nY;
   dimension = probdata->dimension;
   ub = probdata->ub;
   lb = probdata->lb;

   maxnvars = MAX3(4, dimension + 1, ub);

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   SCIP_CALL( SCIPallocBufferArray(scip, &vars, maxnvars) );
   SCIP_CALL( SCIPallocBufferArray(scip, &coeffs, maxnvars) );

   /* create lower bound on objectiver*/
   for (i = 0; i < ub; ++i)
   {
      vars[i] = probdata->isusedvars[i];
      coeffs[i] = 1.0;
   }
   SCIP_CALL( SCIPcreateConsBasicLinear(scip, &(probdata->objbound), "objbound",
         ub, vars, coeffs, (SCIP_Real) lb, SCIPinfinity(scip)) );
   SCIP_CALL( SCIPaddCons(scip, probdata->objbound) );

   /* create constraints ensuring that each inequality is valid for X */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->validconss), ub) );
   for (i = 0; i < ub; ++i)
   {
      for (j = 0; j < dimension; ++j)
         vars[j] = probdata->lhsvars[i][j];
      vars[dimension] = probdata->rhsvars[i];
      coeffs[dimension] = -1.0;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->validconss[i]), nX) );

      for (j = 0; j < nX; ++j)
      {
         int k;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "validcons_%d_%d", i, j);

         for (k = 0; k < dimension; ++k)
            coeffs[k] = probdata->X->points[j][k];

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->validconss[i][j], name, dimension + 1, vars, coeffs, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCons(scip, probdata->validconss[i][j]) );
      }
   }

   /* create constraints ensuring that for each y there is a violated inequality */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->violationconss), nY) );
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "violcons_%d", i);

      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &probdata->violationconss[i], name, ub, probdata->violatedvars[i]) );
      SCIP_CALL( SCIPaddCons(scip, probdata->violationconss[i]) );
   }

   /* link violatedvars and lhs/rhs vars: violatedvars[y][i] = 1 => - lhs * y + rhs <= -eps */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->linkviolconss), nY) );
   coeffs[dimension] = 1.0;
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->linkviolconss[i]), ub) );
      for (j = 0; j < ub; ++j)
      {
         int k;

         for (k = 0; k < dimension; ++k)
         {
            vars[k] = probdata->lhsvars[j][k];
            coeffs[k] = -probdata->Y->points[i][k];
         }
         vars[dimension] = probdata->rhsvars[j];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linkviolcons_%d_%d", i, j);

         SCIP_CALL( SCIPcreateConsBasicIndicator(scip, &probdata->linkviolconss[i][j], name,
               probdata->violatedvars[i][j], dimension + 1, vars, coeffs, -eps) );
         SCIP_CALL( SCIPaddCons(scip, probdata->linkviolconss[i][j]) );
      }
   }

   /* link isusedvars and violatedvars */
   coeffs[0] = 1.0;
   coeffs[1] = -1.0;
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->linkusedconss), nY) );
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->linkusedconss[i]), ub) );
      for (j = 0; j < ub; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "linkusedcons_%d_%d", i, j);

         vars[0] = probdata->violatedvars[i][j];
         vars[1] = probdata->isusedvars[j];

         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->linkusedconss[i][j], name,
               2, vars, coeffs, -1.0, 0.0) );
         SCIP_CALL( SCIPaddCons(scip, probdata->linkusedconss[i][j]) );
      }
   }

   /* add variable upper/lower bound constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->uppervbdconss), ub) );
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->uppervbdconss[i]), dimension) );
      for (j = 0; j < dimension; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "uppervbd_%d_%d", i, j);

         SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &probdata->uppervbdconss[i][j], name,
            probdata->lhsvars[i][j], probdata->isusedvars[i], -1.0, -SCIPinfinity(scip), 0.0) );
         SCIP_CALL( SCIPaddCons(scip, probdata->uppervbdconss[i][j]) );
      }
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->lowervbdconss), ub) );
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->lowervbdconss[i]), dimension) );
      for (j = 0; j < dimension; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lowervbd_%d_%d", i, j);

         SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &probdata->lowervbdconss[i][j], name,
               probdata->lhsvars[i][j], probdata->isusedvars[i], 1.0, 0.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, probdata->lowervbdconss[i][j]) );
      }
   }
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->rhsvbdconss), ub) );
   coeffs[0] = dimension * probdata->absmaxX;
   for (i = 0; i < ub; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "rhsvbd_%d", i);

      SCIP_CALL( SCIPcreateConsBasicVarbound(scip, &probdata->rhsvbdconss[i], name,
            probdata->rhsvars[i], probdata->isusedvars[i], 2 * coeffs[0] , coeffs[0], SCIPinfinity(scip)) );
      SCIP_CALL( SCIPaddCons(scip, probdata->rhsvbdconss[i]) );
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      coeffs[0] = -1.0;
      coeffs[1] = 1.0;

      /* enforce that the isusedvars are sorted */
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->sortusedconss), ub - 1) );
      for (i = 0; i < ub - 1; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sort_isused_%d", i);

         vars[0] = probdata->isusedvars[i];
         vars[1] = probdata->isusedvars[i + 1];

         SCIP_CALL( SCIPcreateConsLinear(scip, &probdata->sortusedconss[i], name, 2, vars, coeffs, -1.0, 0.0,
               TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(scip, probdata->sortusedconss[i]) );
      }

      /* enforce that the violatedvars are sorted */
      SCIP_CALL( SCIPgetBoolParam(scip, "rc/useorbitope", &decision) );
      if ( decision )
      {
#if SCIP_VERSION >= 800
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &probdata->orbitopecons, "orbitope", probdata->violatedvars,
               SCIP_ORBITOPETYPE_FULL, nY, ub, FALSE, FALSE, TRUE, FALSE,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
#else
         SCIP_CALL( SCIPcreateConsOrbitope(scip, &probdata->orbitopecons, "orbitope", probdata->violatedvars,
               SCIP_ORBITOPETYPE_FULL, nY, ub, TRUE, FALSE,
               TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
#endif
         SCIP_CALL( SCIPaddCons(scip, probdata->orbitopecons) );

         /* handle symmetries of Y-points */
         SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision) );
         if ( decision )
         {
            SCIP_VAR** permvars;
            int* perm;
            int npermvars = 0;

            /* turn variable matrix into variable array */
            SCIP_CALL( SCIPallocBufferArray(scip, &permvars, nY * ub) );
            SCIP_CALL( SCIPallocBufferArray(scip, &perm, nY * ub) );

            for (i = 0; i < nY; ++i)
            {
               for (j = 0; j < ub; ++j)
                  permvars[npermvars++] = probdata->violatedvars[i][j];
            }
            assert( npermvars == nY * ub );

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->symresackconss), probdata->nperms) );

            for (i = 0; i < probdata->nperms; ++i)
            {
               int k;
               npermvars = 0;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_%d", i);

               /* construct corresponding permutation */
               for (j = 0; j < nY; ++j)
               {
                  for (k = 0; k < ub; ++k)
                     perm[npermvars++] = probdata->perms[i][j] * ub + k;
               }
               assert( npermvars == nY * ub );

               SCIP_CALL( SCIPcreateSymbreakCons(scip, &probdata->symresackconss[i], name, perm, permvars, npermvars,
                     FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCons(scip, probdata->symresackconss[i]) );
            }

            SCIPfreeBufferArray(scip, &perm);
            SCIPfreeBufferArray(scip, &permvars);
         }
      }
      else
      {
         SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->sortdifflhsconss), ub - 1) );
         coeffs[2] = -2.0;
         coeffs[3] = 2.0;
         for (i = 0; i < ub - 1; ++i)
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sort_diff_%d", i);

            vars[0] = probdata->lhsvars[i][0];
            vars[1] = probdata->lhsvars[i + 1][0];
            vars[2] = probdata->isusedvars[i];
            vars[3] = probdata->isusedvars[i + 1];

            SCIP_CALL( SCIPcreateConsLinear(scip, &probdata->sortdifflhsconss[i], name, 4, vars, coeffs, -SCIPinfinity(scip), 0.0,
                  TRUE, TRUE, TRUE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, probdata->sortdifflhsconss[i]) );
         }
      }
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/hidingsetcuts", &decision) );
   if ( decision )
   {
      int* hidingsetidx;
      int nhidingsetidx;
      int maxnhidingsetidx;
      int cnt;
      int cutcnt = 0;

      SCIP_CALL( computeHidingSetCuts(scip, probdata->X, probdata->Y,
            &hidingsetidx, &nhidingsetidx, &maxnhidingsetidx) );

      assert( nhidingsetidx % 2 == 0 );
      probdata->nhidingsetcuts = ub * nhidingsetidx / 2;

      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->hidingsetcuts, probdata->nhidingsetcuts) );

      for (i = 0; i < ub; ++i)
      {
         cnt = 0;

         for (j = 0; j < nhidingsetidx / 2; ++j)
         {
            (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "hiding_set_%d_%d_inequality_%d",
               hidingsetidx[cnt], hidingsetidx[cnt+1], i);

            vars[0] = probdata->violatedvars[hidingsetidx[cnt++]][i];
            vars[1] = probdata->violatedvars[hidingsetidx[cnt++]][i];

            SCIP_CALL( SCIPcreateConsSetpack(scip, &probdata->hidingsetcuts[cutcnt], name, 2, vars,
                  FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
            SCIP_CALL( SCIPaddCons(scip, probdata->hidingsetcuts[cutcnt++]) );
         }
      }
      assert( cutcnt == probdata->nhidingsetcuts );

      SCIPfreeBlockMemoryArrayNull(scip, &hidingsetidx, maxnhidingsetidx);
   }

   SCIPfreeBufferArray(scip, &coeffs);
   SCIPfreeBufferArray(scip, &vars);

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigRCcompact)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransRCcompact)
{
   int nX;
   int nY;
   int ub;
   int lb;
   int dimension;
   int i;
   SCIP_Bool decision;

   nX = sourcedata->nX;
   nY = sourcedata->nY;
   ub = sourcedata->ub;
   lb = sourcedata->lb;
   dimension = sourcedata->dimension;

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->X, sourcedata->Y, ub, lb, sourcedata->absmaxX,
         sourcedata->perms, sourcedata->nperms, sourcedata->nmaxperms,
         sourcedata->lhsvars, sourcedata->rhsvars, sourcedata->violatedvars, sourcedata->isusedvars,
         sourcedata->validconss, sourcedata->violationconss, sourcedata->linkviolconss, sourcedata->linkusedconss,
         sourcedata->objbound, sourcedata->orbitopecons, sourcedata->symresackconss, sourcedata->sortusedconss,
         sourcedata->sortdifflhsconss, sourcedata->uppervbdconss, sourcedata->lowervbdconss, sourcedata->rhsvbdconss,
         sourcedata->hidingsetcuts, sourcedata->nhidingsetcuts) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->objbound, &(*targetdata)->objbound) );
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, nX, (*targetdata)->validconss[i], (*targetdata)->validconss[i]) );
   }
   SCIP_CALL( SCIPtransformConss(scip, nY, (*targetdata)->violationconss, (*targetdata)->violationconss) );
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, ub, (*targetdata)->linkviolconss[i], (*targetdata)->linkviolconss[i]) );
   }
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, ub, (*targetdata)->linkusedconss[i], (*targetdata)->linkusedconss[i]) );
   }
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, dimension, (*targetdata)->uppervbdconss[i], (*targetdata)->uppervbdconss[i]) );
   }
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, dimension, (*targetdata)->lowervbdconss[i], (*targetdata)->lowervbdconss[i]) );
   }
   SCIP_CALL( SCIPtransformConss(scip, ub, (*targetdata)->rhsvbdconss, (*targetdata)->rhsvbdconss) );

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "rc/useorbitope", &decision) );
      if ( decision )
      {
         SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->orbitopecons, &(*targetdata)->orbitopecons) );

         SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision) );
         if ( decision )
         {
            SCIP_CALL( SCIPtransformConss(scip, (*targetdata)->nperms, (*targetdata)->symresackconss, (*targetdata)->symresackconss) );
         }
      }
      else
      {
         SCIP_CALL( SCIPtransformConss(scip, ub - 1, (*targetdata)->sortdifflhsconss, (*targetdata)->sortdifflhsconss) );
      }
      SCIP_CALL( SCIPtransformConss(scip, ub - 1, (*targetdata)->sortusedconss, (*targetdata)->sortusedconss) );
   }

   SCIP_CALL( SCIPtransformConss(scip, sourcedata->nhidingsetcuts, (*targetdata)->hidingsetcuts, (*targetdata)->hidingsetcuts) );

   /* transform all variables */
   for (i = 0; i < ub; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, dimension, (*targetdata)->lhsvars[i], (*targetdata)->lhsvars[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, ub, (*targetdata)->rhsvars, (*targetdata)->rhsvars) );
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, ub, (*targetdata)->violatedvars[i], (*targetdata)->violatedvars[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, ub, (*targetdata)->isusedvars, (*targetdata)->isusedvars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransRCcompact)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolRCcompact)
{
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolRCcompact)
{  /*lint --e{715}*/
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateCompact(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX,            /**< maximum absolute value of a coordinate in X */
   SCIP_Real**           inequalities,       /**< allocated array to store inequalities (or NULL if not needed) */
   int**                 separatedpoints,    /**< allocated array to store separated points per inequality
                                                (or NULL if not needed) */
   int*                  nseparatedpoints,   /**< allocated array to store number of separated points per inequality
                                                (or NULL if not needed) */
   int                   ninequalities,      /**< number of inequalities encoded in previous data structures
                                                (or -1 if not needed) */
   int                   maxninequalities    /**< maximum number of inequalities that can be stored
                                                (if allocated) */
   )
{
   SCIP_PROBDATA* probdata;
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   SCIP_Bool decision;
   SCIP_Bool decision2;
   SCIP_Bool ubcomputed = FALSE;
   int** perms = NULL;
   int nperms = 0;
   int nmaxperms = 0;
   int i;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub != NULL );
   assert( *ub > 0 );
   assert( absmaxX > 0 );
   assert( inequalities != NULL || ninequalities == -1 );
   assert( separatedpoints != NULL || ninequalities == -1 );
   assert( nseparatedpoints != NULL || ninequalities == -1 );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigRCcompact) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransRCcompact) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransRCcompact) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolRCcompact) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolRCcompact) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );

   /* compute upper bound if trivial upper bound is given */
   if ( *ub == INT_MAX )
   {
      dd_set_global_constants();
      generators = constructGeneratorMatrixPoints(X, NULL, 0);

      facetsconvexhull = computeConvexHullFacets(scip, generators, &ubcomputed);

      dd_FreeMatrix(generators);

      if ( ubcomputed )
         *ub = facetsconvexhull->rowsize;
      else
         dd_free_global_constants();
   }

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision2) );
   if ( decision && decision2 )
   {
      SCIP_CALL( computeSymmetries(scip, X, Y, &perms, &nperms, &nmaxperms) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, X, Y, *ub, lb,  absmaxX, perms, nperms, nmaxperms,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0) );

   SCIP_CALL( SCIPcreateVariables(scip, probdata, absmaxX) );
   SCIP_CALL( SCIPcreateConstraints(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   /* free original permutations (are copied by probdataCreate()) */
   for (i = 0; i < nperms; ++i)
   {
      SCIPfreeBlockMemoryArrayNull(scip, &perms[i], Y->ndatapoints);
   }
   SCIPfreeBlockMemoryArrayNull(scip, &perms, nmaxperms);

   /* compute initial solution */
   if ( ubcomputed || ninequalities > 0 )
   {
      SCIP_SOL* sol;
      SCIP_Bool stored;

      SCIP_CALL( SCIPcreateOrigSol(scip, &sol, NULL) );

      if ( ubcomputed )
      {
         SCIP_CALL( SCIPgetSolutionCompactModelCDD(scip, probdata, sol, facetsconvexhull) );
      }
      else
      {
         SCIP_CALL( SCIPgetSolutionCompactModelExplicit(scip, probdata, sol, inequalities,
               separatedpoints, nseparatedpoints, ninequalities) );
      }

      SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      dd_free_global_constants();
   }

   return SCIP_OKAY;
}

/** return variables encoding right-hand sides of relaxation */
SCIP_VAR** SCIPprobdataGetRhsvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->rhsvars;
}

/** return variables encoding left-hand sides of relaxation */
SCIP_VAR***  SCIPprobdataGetLhsvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->lhsvars;
}

/** return variables indicating whether an inequality is used */
SCIP_VAR**  SCIPprobdataGetUsedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->isusedvars;
}

/** return variables indicating whether an inequality is violated */
SCIP_VAR***  SCIPprobdataGetViolatedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->violatedvars;
}

/** return constraints linking violation of variables and inequality */
SCIP_CONS***  SCIPprobdataGetLinkviolconss(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->linkviolconss;
}

/** return the set of feasible points */
Datapoints*  SCIPprobdataCompactGetX(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->X;
}

/** return the set of infeasible points */
Datapoints*  SCIPprobdataCompactGetY(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->Y;
}

/** return the upper bound on needed inequalities */
int SCIPprobdataGetUb(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );

   return probdata->ub;
}

/**@} */
