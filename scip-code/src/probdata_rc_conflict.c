/**@file   probdata_rc_conflict.c
 * @brief  Problem data for computing RC using a conflict based model
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "compute_symmetries.h"
#include "cons_conflict.h"
#include "hiding_sets.h"
#include "datapoints.h"
#include "probdata_rc_conflict.h"
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
#include "scip/cons_orbitope.h"
#include "scip/cons_setppc.h"
#include "scip/cons_symresack.h"

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
   SCIP_VAR***           violatedvars;       /**< variables indicating whether constraint is violated */
   SCIP_VAR**            isusedvars;         /**< variables indicating whether constraint is used */
   SCIP_CONS**           violationconss;     /**< constraints ensuring each infeasible point violates a constraint */
   SCIP_CONS***          linkusedconss;      /**< constraints to link isusedvars with constraints */
   SCIP_CONS*            objbound;           /**< constraint modeling lower bound on objective */
   SCIP_CONS*            conflictcons;       /**< constraint modeling the conflicts of the model */

   SCIP_CONS*            orbitopecons;       /**< constraint to handle permutations of inequality classes */
   SCIP_CONS**           symresackconss;     /**< constraints to handle permutations of Y */
   SCIP_CONS**           sortconss;          /**< constraints sorting isusedvars */
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
   SCIP_VAR***           violatedvars,       /**< variables indicating whether constraint is violated (Y x ub)*/
   SCIP_VAR**            isusedvars,         /**< variables indicating whether constraint is used */
   SCIP_CONS**           violationconss,     /**< constraints ensuring each infeasible point violated a constraint */
   SCIP_CONS***          linkusedconss,      /**< constraints to link isusedvars with constraints (Y x ub) */
   SCIP_CONS*            objbound,           /**< constraint modeling lower bound on objective */
   SCIP_CONS*            conflictcons,       /**< constraint modeling conflicts of the model */
   SCIP_CONS*            orbitopecons,       /**< constraint to handle permutations of inequality classes */
   SCIP_CONS**           symresackconss,     /**< constraints to handle permutations of Y */
   SCIP_CONS**           sortconss,          /**< constraints sorting isusedvars */
   SCIP_CONS**           hidingsetcuts,      /**< hiding set based cuts */
   int                   nhidingsetcuts      /**< number of hiding set cuts stored in hidingsetcuts */
   )
{
   int nX;
   int nY;
   int dimension;
   int i;

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
   (*probdata)->conflictcons = conflictcons;
   (*probdata)->orbitopecons = orbitopecons;

   /* possible copy variable arrays */
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
   if ( violationconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->violationconss, violationconss, nY) );
   }
   else
      (*probdata)->violationconss = NULL;

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

   if ( sortconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->sortconss, sortconss, ub - 1) );
   }
   else
      sortconss = NULL;
   if ( symresackconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->symresackconss, symresackconss, nperms) );
   }
   else
      (*probdata)->symresackconss = NULL;

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
   int nY;
   int ub;
   SCIP_Bool decision;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->X != NULL );
   assert( (*probdata)->Y != NULL );
   assert( (*probdata)->ub > 0 );

   nY = (*probdata)->nY;
   ub = (*probdata)->ub;

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

   /* release all constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->conflictcons) );
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

   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->violationconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->violationconss, nY);

   for (i = 0; i < (*probdata)->nhidingsetcuts; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->hidingsetcuts[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->hidingsetcuts, (*probdata)->nhidingsetcuts);

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      SCIP_CALL( SCIPgetBoolParam(scip, "rc/useorbitope", &decision) );

      if ( decision )
      {
         SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->orbitopecons) );
         for (i = 0; i < ub - 1; ++i)
         {
            SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->sortconss[i]) );
         }
         SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->sortconss, ub - 1);

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
   }

   /* free memory of arrays */
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->isusedvars, ub);
   for (i = nY - 1; i >= 0; --i)
   {
      SCIPfreeBlockMemoryArray(scip, &(*probdata)->violatedvars[i], ub);
   }
   SCIPfreeBlockMemoryArray(scip, &(*probdata)->violatedvars, nY);
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
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   int nY;
   int ub;
   int i;
   int j;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( probdata->nY > 0 );
   assert( probdata->ub > 0 );

   nY = probdata->nY;
   ub = probdata->ub;

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

      SCIP_CALL( SCIPcreateVarBasic(scip, &probdata->isusedvars[i], name, 0.0, 1.0, 1.0, SCIP_VARTYPE_BINARY) );
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
   int nY;
   int ub;
   int lb;
   int dimension;
   int i;
   int j;
   SCIP_VAR** vars;
   SCIP_Real* coeffs;
   int maxnvars;
   SCIP_Bool decision;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->Y != NULL );
   assert( probdata->Y->points != NULL );
   assert( probdata->nY > 0 );
   assert( probdata->ub > 0 );

   nY = probdata->nY;
   ub = probdata->ub;
   lb = probdata->lb;
   dimension = probdata->dimension;

   maxnvars = MAX(2, ub);

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

   /* create constraints ensuring that for each y there is exactly one violated inequality */
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->violationconss), nY) );
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "violcons_%d", i);

      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &probdata->violationconss[i], name, ub, probdata->violatedvars[i]) );
      SCIP_CALL( SCIPaddCons(scip, probdata->violationconss[i]) );
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

   SCIP_CALL( SCIPcreateConsBasicConflict(scip, &probdata->conflictcons, "conflict",
         probdata->X, probdata->Y, ub, probdata->violatedvars, probdata->absmaxX) );
   SCIP_CALL( SCIPaddCons(scip, probdata->conflictcons) );

   /* possibly handle symmetries */
   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &probdata->sortconss, ub - 1) );
      for (i = 0; i < ub - 1; ++i)
      {
         vars[0] = probdata->isusedvars[i];
         vars[1] = probdata->isusedvars[i + 1];

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sort_%d", i);
         SCIP_CALL( SCIPcreateConsBasicLinear(scip, &probdata->sortconss[i], name,
               2, vars, coeffs, 0.0, SCIPinfinity(scip)) );
         SCIP_CALL( SCIPaddCons(scip, probdata->sortconss[i]) );
      }

      SCIP_CALL( SCIPgetBoolParam(scip, "rc/useorbitope", &decision) );
      if ( decision )
      {
#if SCIP_VERSION >= 800
         SCIP_CALL( SCIPcreateConsBasicOrbitope(scip, &probdata->orbitopecons, "orbitope",
               probdata->violatedvars, SCIP_ORBITOPETYPE_FULL, nY, ub, FALSE, FALSE, TRUE, FALSE) );
#else
	 SCIP_CALL( SCIPcreateConsBasicOrbitope(scip, &probdata->orbitopecons, "orbitope",
               probdata->violatedvars, SCIP_ORBITOPETYPE_FULL, nY, ub, TRUE, FALSE) );
#endif
         SCIP_CALL( SCIPaddCons(scip, probdata->orbitopecons) );

         SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision) );
         if ( decision )
         {
            SCIP_VAR** permvars;
            int* perm;
            int npermvars = 0;

            /* turn variable matrix into variable array */
            SCIP_CALL( SCIPallocBufferArray(scip, &permvars, nY * dimension) );
            SCIP_CALL( SCIPallocBufferArray(scip, &perm, nY * dimension) );

            for (i = 0; i < nY; ++i)
            {
               for (j = 0; j < dimension; ++j)
                  permvars[npermvars++] = probdata->violatedvars[i][j];
            }
            assert( npermvars == nY * dimension );

            SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(probdata->symresackconss), probdata->nperms) );

            for (i = 0; i < probdata->nperms; ++i)
            {
               int k;
               npermvars = 0;

               (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "symresack_%d", i);

               /* construct corresponding permutation */
               for (j = 0; j < nY; ++j)
               {
                  for (k = 0; k < dimension; ++k)
                     perm[npermvars++] = probdata->perms[i][j] * dimension + k;
               }
               assert( npermvars == nY * dimension );

               SCIP_CALL( SCIPcreateSymbreakCons(scip, &probdata->symresackconss[i], name, perm, permvars, npermvars,
                     FALSE, TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
               SCIP_CALL( SCIPaddCons(scip, probdata->symresackconss[i]) );
            }

            SCIPfreeBufferArray(scip, &perm);
            SCIPfreeBufferArray(scip, &permvars);
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
SCIP_DECL_PROBDELORIG(probdelorigRCconflict)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransRCconflict)
{
   int nY;
   int ub;
   int lb;
   int i;
   SCIP_Bool decision;

   nY = sourcedata->nY;
   ub = sourcedata->ub;
   lb = sourcedata->lb;

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->X, sourcedata->Y, ub, lb, sourcedata->absmaxX,
         sourcedata->perms, sourcedata->nperms, sourcedata->nmaxperms,
         sourcedata->violatedvars, sourcedata->isusedvars, sourcedata->violationconss, sourcedata->linkusedconss,
         sourcedata->objbound, sourcedata->conflictcons, sourcedata->orbitopecons, sourcedata->symresackconss,
         sourcedata->sortconss, sourcedata->hidingsetcuts, sourcedata->nhidingsetcuts) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->objbound, &(*targetdata)->objbound) );
   SCIP_CALL( SCIPtransformCons(scip, (*targetdata)->conflictcons, &(*targetdata)->conflictcons) );
   SCIP_CALL( SCIPtransformConss(scip, nY, (*targetdata)->violationconss, (*targetdata)->violationconss) );
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPtransformConss(scip, ub, (*targetdata)->linkusedconss[i], (*targetdata)->linkusedconss[i]) );
   }
   SCIP_CALL( SCIPtransformConss(scip, sourcedata->nhidingsetcuts, (*targetdata)->hidingsetcuts, (*targetdata)->hidingsetcuts) );

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/handlesymmetry", &decision) );
   if ( decision )
   {
      SCIP_CALL( SCIPtransformConss(scip, ub - 1, (*targetdata)->sortconss, (*targetdata)->sortconss) );

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
   }

   /* transform all variables */
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPtransformVars(scip, ub, (*targetdata)->violatedvars[i], (*targetdata)->violatedvars[i]) );
   }
   SCIP_CALL( SCIPtransformVars(scip, ub, (*targetdata)->isusedvars, (*targetdata)->isusedvars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransRCconflict)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolRCconflict)
{
   assert(probdata != NULL);

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolRCconflict)
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
SCIP_RETCODE SCIPprobdataCreateConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX             /**< maximum absolute value of a coordinate in X */
   )
{
   SCIP_PROBDATA* probdata;
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   int** perms = NULL;
   int nperms = 0;
   int nmaxperms = 0;
   int i;
   SCIP_Bool decision;
   SCIP_Bool ubcomputed = FALSE;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub != NULL );
   assert( *ub > 0 );
   assert( absmaxX > 0 );

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigRCconflict) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransRCconflict) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransRCconflict) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolRCconflict) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolRCconflict) );

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

   SCIP_CALL( SCIPgetBoolParam(scip, "rc/usesymresacks", &decision) );
   if ( decision )
   {
      SCIP_CALL( computeSymmetries(scip, X, Y, &perms, &nperms, &nmaxperms) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, X, Y, *ub, lb, absmaxX, perms, nperms, nmaxperms,
         NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0) );

   SCIP_CALL( SCIPcreateVariables(scip, probdata) );
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
   if ( ubcomputed )
   {
      SCIP_SOL* sol;
      SCIP_Bool stored;

      SCIP_CALL( SCIPcreateOrigSol(scip, &sol, NULL) );

      SCIP_CALL( SCIPgetSolutionConflictModel(scip, probdata, sol, facetsconvexhull) );

      SCIP_CALL( SCIPaddSol(scip, sol, &stored) );
      SCIP_CALL( SCIPfreeSol(scip, &sol) );

      dd_free_global_constants();
   }

   return SCIP_OKAY;
}

/**< returns feasible points */
Datapoints* SCIPprobdataConflictGetX(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );
   return probdata->X;
}

/**< returns infeasible points */
Datapoints* SCIPprobdataConflictGetY(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );
   return probdata->Y;
}

/**< returns upper bound on number of classes */
int SCIPprobdataConflictGetUb(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );
   return probdata->ub;
}

/**< returns violation variables */
SCIP_VAR*** SCIPprobdataConflictGetViolatedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );
   return probdata->violatedvars;
}

/**< returns isused variables */
SCIP_VAR** SCIPprobdataConflictGetIsusedvars(
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   assert( probdata != NULL );
   return probdata->isusedvars;
}

/**@} */
