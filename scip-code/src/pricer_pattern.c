/*
  This file is based on the file pricer_binpack.c distributed by via SCIP:

  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
  *                                                                           *
  *                  This file is part of the program and library             *
  *         SCIP --- Solving Constraint Integer Programs                      *
  *                                                                           *
  *    Copyright (C) 2002-2021 Zuse Institute Berlin                          *
  *                                                                           *
  *                                                                           *
  *  SCIP is distributed under the terms of the ZIB Academic License.         *
  *                                                                           *
  *  You should have received a copy of the ZIB Academic License              *
  *  along with SCIP; see the file COPYING. If not visit scipopt.org.         *
  *                                                                           *
  * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

  @file   pricer_binpacking.c
  @brief  Binpacking variable pricer
  @author Timo Berthold
  @author Stefan Heinz

  The source code is adapted to our problem data structures and creates the pricing
  problem in a slightly different way. Except for this, the source code remains
  unchanged.
*/



/**@file   pricer_pattern.c
 * @brief  stable set variable pricer
 * @author Timo Berthold
 * @author Stefan Heinz
 * @author Christopher Hojny
 *
 * Implementation of the variable pricer. To detect a new pattern, we
 * search for a maximum weight (w.r.t. dual values) set of infeasible points Y
 * that can be separated by a single linear inequality from X. If we are in
 * a branching node, we respect the branching decisions encoded by SAME/DIFFER
 * constraints. Moreover, if a variable is fixed to 0, this decision is also
 * incorporated into the pricing problem to avoid generating this variable again.
 * To implement these latter two changes, we use slight modifications of the methods
 * addBranchingDecisionConss() and addFixedVarsConss() of the Binpacking example
 * in SCIP, which are defined in pricer_binpacking.c. For this reason, the authors
 * of that file are also added as authors of this file.
 */

#include "cons_samediff.h"
#include "datapoints.h"
#include "hiding_sets.h"
#include "maximal_separation.h"
#include "pricer_pattern.h"
#include "probdata_rc_cg.h"
#include "vardata_binpacking.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"
#include "scip/cons_setppc.h"

#include <assert.h>

#define PRICER_NAME            "pattern"
#define PRICER_DESC            "pattern variable pricer"
#define PRICER_PRIORITY        0
#define PRICER_DELAY           TRUE     /* only call pricer if all problem variables have non-negative reduced costs */


struct SCIP_PricerData
{
   SCIP_CONS**           coverconss;         /**< covering constraints of main problem */
   int                   ncoverconss;        /**< number of covering constraint */
   Datapoints*           X;                  /**< set of feasible points */
   Datapoints*           Y;                  /**< set of infeasible points */
   int                   nX;                 /**< number of feasible points */
   int                   nY;                 /**< number of infeasible points */
   int                   dimension;          /**< dimension of points in X and Y */
   int                   absmaxX;            /**< maximum absolute value of coordinate in X */
   SCIP_CONSHDLR*        conshdlr;           /**< pointer to same/diff conshldr */

   SCIP*                 subscip;            /**< SCIP instance of pricing problem */
   SCIP_VAR**            selectvars;         /**< variable for selecting y's in pricing problem */
   SCIP_VAR**            negselectvars;      /**< negation of variable for selecting y's in pricing problem */
   SCIP_VAR**            ineqvars;           /**< variables encoding an inequality in type (a,-b) <= 0 */
   SCIP_CONS**           negconss;           /**< constraints defining negated variables */
   SCIP_CONS**           validconss;         /**< constraints for finding valid inequality in pricing problem */
   SCIP_CONS**           sepaconss;          /**< constraints for finding separating inequality in pricing problem */
   SCIP_CONS**           negsepaconss;       /**< constraints for forbidding separating inequality in pricing problem */
   SCIP_CONS**           hidingsetcuts;      /**< hiding set cuts to strengthen pricing problem */
   int*                  hidingsetidx;       /**< array containing pairs of infeasible points encoding hiding set */
   int                   nhidingsetidx;      /**< number of indices stores in hidingsetidx */
   int                   maxnhidingsetidx;   /**< maximum number of indices that can be stored in hidingsetidx */
   int                   nhidingsetcuts;     /**< number of hiding set cuts */
};

/*
 * Local methods
 */

/** add branching decisions constraints to the sub SCIP */
static
SCIP_RETCODE addBranchingDecisionConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array of the subscuip oder variables */
   SCIP_CONSHDLR*        conshdlr            /**< constraint handler for branching data */
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


/** avoid to generate columns which are fixed to zero; therefore add for each variable which is fixed to zero a
 *  corresponding logicor constraint to forbid this column
 */
static
SCIP_RETCODE addFixedVarsConss(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP*                 subscip,            /**< pricing SCIP data structure */
   SCIP_VAR**            vars,               /**< variable array of the subscuip */
   SCIP_CONS**           conss,              /**< array of setppc constraint for each item one */
   int                   nitems              /**< number of items */
   )
{
   SCIP_VAR** origvars;
   int norigvars;

   SCIP_CONS* cons;
   int* consids;
   int nconsids;
   int consid;
   int nvars;

   SCIP_VAR** logicorvars;
   SCIP_VAR* var;
   SCIP_VARDATA* vardata;
   SCIP_Bool needed;
   int nlogicorvars;

   int v;
   int c;
   int o;

   /* collect all variable which are currently existing */
   origvars = SCIPgetVars(scip);
   norigvars = SCIPgetNVars(scip);

   /* loop over all these variables and check if they are fixed to zero */
   for (v = 0; v < norigvars; ++v)
   {
      assert(SCIPvarGetType(origvars[v]) == SCIP_VARTYPE_BINARY);

      /* if the upper bound is smaller than 0.5 it follows due to the integrality that the binary variable is fixed to zero */
      if ( SCIPvarGetUbLocal(origvars[v]) < 0.5 )
      {
         /* collect the constraints/patterns the variable belongs to */
         vardata = SCIPvarGetData(origvars[v]);
         nconsids = SCIPvardataGetNConsids(vardata);
         consids = SCIPvardataGetConsids(vardata);
         needed = TRUE;

         SCIP_CALL( SCIPallocBufferArray(subscip, &logicorvars, nitems) );
         nlogicorvars = 0;
         consid = consids[0];
         nvars = 0;

         /* loop over these items and create a linear (logicor) constraint which forbids this item combination in the
          * pricing problem; thereby check if this item combination is already forbidden
          */
         for (c = 0, o = 0; o < nitems && needed; ++o)
         {
            assert( o <= consid );
            cons = conss[o];

            /* if the constraint is already disabled */
            if ( o == consid && !SCIPconsIsEnabled(cons) )
            {
               needed = FALSE;
               break;
            }

            var = vars[nvars];
            nvars++;
            assert(var != NULL);

            if ( o == consid )
            {
               SCIP_CALL( SCIPgetNegatedVar(subscip, var, &var) );
            }

            logicorvars[nlogicorvars] = var;
            nlogicorvars++;

            if ( o == consid )
            {
               c++;
               if ( c == nconsids )
                  consid = nitems + 100;
               else
               {
                  assert(consid < consids[c]);
                  consid = consids[c];
               }
            }
         }

         if ( needed )
         {
            SCIP_CALL( SCIPcreateConsBasicLogicor(subscip, &cons, SCIPvarGetName(origvars[v]), nlogicorvars, logicorvars) );
            SCIP_CALL( SCIPsetConsInitial(subscip, cons, FALSE) );

            SCIP_CALL( SCIPaddCons(subscip, cons) );
            SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
         }

         SCIPfreeBufferArray(subscip, &logicorvars);
      }
   }

   return SCIP_OKAY;
}


/** creates model for pricing new variables */
static
SCIP_RETCODE createPricingModel(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP*                 subscip,            /**< SCIP instance of the pricing model */
   SCIP_PRICERDATA*      pricerdata,         /**< pointer to pricer data */
   SCIP_Bool             atrootnode          /**< whether we create the pricing problem at the root node */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VAR** tmpvars;
   SCIP_Real* tmpvals;
   SCIP_Real** Xpoints;
   SCIP_Real** Ypoints;
   SCIP_Real eps;
   int nX;
   int nY;
   int absmaxX;
   int dimension;
   int i;
   int j;
   SCIP_Bool decision;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( pricerdata != NULL );
   assert( pricerdata->X != NULL );
   assert( pricerdata->Y != NULL );
   assert( pricerdata->X->points != NULL );
   assert( pricerdata->Y->points != NULL );
   assert( pricerdata->absmaxX > 0 );

   Xpoints = pricerdata->X->points;
   Ypoints = pricerdata->Y->points;
   nX = pricerdata->nX;
   nY = pricerdata->nY;
   dimension = pricerdata->dimension;
   absmaxX = pricerdata->absmaxX;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* basic setup */
   SCIP_CALL( SCIPcreateProbBasic(subscip, "setGenerator") );
   if ( ! atrootnode )
   {
      SCIP_CALL( SCIPsetEmphasis(subscip, SCIP_PARAMEMPHASIS_FEASIBILITY, TRUE) );
   }
   SCIP_CALL( SCIPsetObjsense(subscip, SCIP_OBJSENSE_MAXIMIZE) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* load plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->selectvars, nY) );
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "select_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &pricerdata->selectvars[i], name, 0.0, 1.0, 0.0,
            SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, pricerdata->selectvars[i]) );
   }

   if ( ! atrootnode )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->negselectvars, nY) );
      for (i = 0; i < nY; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "negselect_%d", i);

         SCIP_CALL( SCIPcreateVar(subscip, &pricerdata->negselectvars[i], name, 0.0, 1.0, 0.0,
               SCIP_VARTYPE_BINARY, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
         SCIP_CALL( SCIPaddVar(subscip, pricerdata->negselectvars[i]) );
      }
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->ineqvars, dimension + 1) );
   for (i = 0; i < dimension; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "a_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &pricerdata->ineqvars[i], name, -1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, pricerdata->ineqvars[i]) );
   }
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b");

   SCIP_CALL( SCIPcreateVar(subscip, &pricerdata->ineqvars[dimension], name,
         -dimension * absmaxX, dimension * absmaxX, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, pricerdata->ineqvars[dimension]) );

   /* create constraints */
   SCIP_CALL( SCIPallocBufferArray(scip, &tmpvals, dimension + 1 ) );

   if ( ! atrootnode )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->negconss, nY) );
      SCIP_CALL( SCIPallocBufferArray(subscip, &tmpvars, 2) );
      for (i = 0; i < nY; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "negatecons_%d", i);

         tmpvars[0] = pricerdata->selectvars[i];
         tmpvars[1] = pricerdata->negselectvars[i];

         SCIP_CALL( SCIPcreateConsBasicSetpart(subscip, &pricerdata->negconss[i], name, 2, tmpvars) );
         SCIP_CALL( SCIPaddCons(subscip, pricerdata->negconss[i]) );
      }
      SCIPfreeBufferArray(subscip, &tmpvars);
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->validconss, nX) );
   tmpvals[dimension] = -1.0;
   for (i = 0; i < nX; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "validcons_%d", i);

      for (j = 0; j < dimension; ++j)
         tmpvals[j] = Xpoints[i][j];
      tmpvals[dimension] = -1.0;

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &pricerdata->validconss[i], name,
            dimension + 1, pricerdata->ineqvars, tmpvals, -SCIPinfinity(subscip), 0.0) );
      SCIP_CALL( SCIPaddCons(subscip, pricerdata->validconss[i]) );
   }

   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->sepaconss, nY) );
   tmpvals[dimension] = 1.0;
   for (i = 0; i < nY; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sepacons_%d", i);

      for (j = 0; j < dimension; ++j)
         tmpvals[j] = -Ypoints[i][j];

      SCIP_CALL( SCIPcreateConsBasicIndicator(subscip, &pricerdata->sepaconss[i], name,
            pricerdata->selectvars[i], dimension + 1, pricerdata->ineqvars, tmpvals, -eps) );
      SCIP_CALL( SCIPaddCons(subscip, pricerdata->sepaconss[i]) );
   }

   if ( ! atrootnode )
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->negsepaconss, nY) );
      tmpvals[dimension] = -1.0;
      for (i = 0; i < nY; ++i)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "negsepacons_%d", i);

         for (j = 0; j < dimension; ++j)
            tmpvals[j] = Ypoints[i][j];

         SCIP_CALL( SCIPcreateConsBasicIndicator(subscip, &pricerdata->negsepaconss[i], name,
               pricerdata->negselectvars[i], dimension + 1, pricerdata->ineqvars, tmpvals, 0.0) );
         SCIP_CALL( SCIPaddCons(subscip, pricerdata->negsepaconss[i]) );
      }
   }

   /* check whether we have to add hiding set cuts */
   SCIP_CALL( SCIPgetBoolParam(scip, "rc/hidingsetcuts", &decision) );
   if ( decision )
   {
      SCIP_VAR* vars[2];
      int* hidingsetidx;
      int cnt;

      if ( pricerdata->nhidingsetidx == 0 )
      {
         SCIP_CALL( computeHidingSetCuts(subscip, pricerdata->X, pricerdata->Y,
               &pricerdata->hidingsetidx, &pricerdata->nhidingsetidx, &pricerdata->maxnhidingsetidx) );
      }

      assert( pricerdata->nhidingsetidx % 2 == 0 );
      pricerdata->nhidingsetcuts = pricerdata->nhidingsetidx / 2;

      SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &pricerdata->hidingsetcuts, pricerdata->nhidingsetcuts) );

      hidingsetidx = pricerdata->hidingsetidx;
      cnt = 0;
      for (j = 0; j < pricerdata->nhidingsetcuts; ++j)
      {
         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "hiding_set_%d_%d", hidingsetidx[cnt], hidingsetidx[cnt+1]);

         vars[0] = pricerdata->selectvars[hidingsetidx[cnt++]];
         vars[1] = pricerdata->selectvars[hidingsetidx[cnt++]];

         SCIP_CALL( SCIPcreateConsSetpack(subscip, &pricerdata->hidingsetcuts[j], name, 2, vars,
               FALSE, TRUE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );
         SCIP_CALL( SCIPaddCons(subscip, pricerdata->hidingsetcuts[j]) );
      }
   }

   SCIPfreeBufferArray(scip, &tmpvals);

   return SCIP_OKAY;
}


/** frees data of pricing model */
static
SCIP_RETCODE freePricingModel(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP*                 subscip,            /**< SCIP instance of the pricing model */
   SCIP_PRICERDATA*      pricerdata,         /**< pointer to pricer data */
   SCIP_Bool             atrootnode          /**< whether we free the pricing problem at the root node */
   )
{
   int nX;
   int nY;
   int dimension;
   int i;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( pricerdata != NULL );
   assert( pricerdata->X != NULL );
   assert( pricerdata->Y != NULL );

   nX = pricerdata->nX;
   nY = pricerdata->nY;
   dimension = pricerdata->dimension;

   /* release variables */
   if ( ! atrootnode )
   {
      for (i = 0; i < nY ; ++i)
      {
         SCIP_CALL( SCIPreleaseVar(subscip, &pricerdata->negselectvars[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->negselectvars, nY);
   }
   for (i = 0; i < nY ; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &pricerdata->selectvars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->selectvars, nY);
   for (i = 0; i < dimension + 1; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &pricerdata->ineqvars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->ineqvars, dimension + 1);

   /* release conss */
   if ( ! atrootnode )
   {
      for (i = 0; i < nY; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(subscip, &pricerdata->negsepaconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->negsepaconss, nY);
      for (i = 0; i < nY; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(subscip, &pricerdata->negconss[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->negconss, nY);
   }
   for (i = 0; i < nX; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &pricerdata->validconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->validconss, nX);
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &pricerdata->sepaconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->sepaconss, nY);

   if ( pricerdata->nhidingsetcuts > 0 )
   {
      for (i = 0; i < pricerdata->nhidingsetcuts; ++i)
      {
         SCIP_CALL( SCIPreleaseCons(subscip, &pricerdata->hidingsetcuts[i]) );
      }
      SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->hidingsetcuts, pricerdata->nhidingsetcuts);
   }

   /* free memory */
   if ( pricerdata->maxnhidingsetidx != 0 )
   {
      SCIPfreeBlockMemoryArrayNull(subscip, &pricerdata->hidingsetidx, pricerdata->maxnhidingsetidx);
   }
   pricerdata->nhidingsetcuts = 0;
   pricerdata->nhidingsetidx = 0;
   pricerdata->maxnhidingsetidx = 0;

   /* free problem */
   SCIP_CALL( SCIPfreeTransform(subscip) );
   SCIP_CALL( SCIPfree(&subscip) );

   return SCIP_OKAY;
}


/** updates objective of pricing model by current dual solution */
static
SCIP_RETCODE changeObjectivePricingProblem(
   SCIP*                 scip,               /**< main SCIP instance */
   SCIP*                 subscip,            /**< SCIP instance of pricing model */
   SCIP_CONS**           coverconss,         /**< covering constraints */
   int                   ncoverconss,        /**< number of covering constraints */
   SCIP_VAR**            dualvars            /**< relevant variables of pricing model */
   )
{
   int c;
   SCIP_CONS* cons;
   SCIP_Real dual;

   assert( scip != NULL );
   assert( subscip != NULL );
   assert( coverconss != NULL );
   assert( ncoverconss > 0 );

   /* per constraint, get the dual multipliers */
   for (c = 0; c < ncoverconss; ++c)
   {
      cons = coverconss[c];

      /* get dual multiplier */
      dual = SCIPgetDualsolSetppc(scip, cons);

      /* adapt objective in separation problem */
      SCIP_CALL( SCIPchgVarObj(subscip, dualvars[c], dual) );
   }

   return SCIP_OKAY;
}


/** solves the pricing problem */
static
SCIP_RETCODE SCIPsolvePricingProblem(
   SCIP*                 scip,               /**< main SCIP instance */
   SCIP*                 subscip,            /**< SCIP instance of pricing problem */
   Datapoints*           X,                  /**< set of feasible points */
   int                   nX,                 /**< number of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int                   nY,                 /**< number of infeasible points */
   int                   dimension,          /**< ambient dimension of infeasible points */
   int                   absmaxX,            /**< maximum absolute value of coordinate in X */
   SCIP_VAR**            ineqvars,           /**< array of variables encoding inequality defining priced pattern */
   SCIP_CONS**           coverconss,         /**< covering constraints of main problem */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler for branching data */
   SCIP_Bool*            addvar              /**< pointer to store whether a new variable coudl be added */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_VARDATA* vardata;
   SCIP_PROBDATA* probdata;
   SCIP_SOL* subscipsol;
   SCIP_VAR* newvar;
   SCIP_Real* separatedineq;
   SCIP_Real subscipobj;
   SCIP_Real val;
   SCIP_Real eps;
   SCIP_Bool maxsepa;
   int* consids;
   int nconss;
   int i;
   int j;

   SCIP_Real cutoffbound;
   int sollimit;

   *addvar = FALSE;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/cutoffepsilon", &cutoffbound) );
   SCIP_CALL( SCIPgetIntParam(scip, "rc/maxsolpricer", &sollimit) );
   SCIP_CALL( SCIPgetBoolParam(scip, "rc/maxsepa", &maxsepa) );
   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* solve sub SCIP to find new patterns/sets for master problem  */
   SCIP_CALL( SCIPsetObjlimit(subscip, 1.0 + cutoffbound) );
   SCIP_CALL( SCIPsetIntParam(subscip, "limits/bestsol", sollimit) );

   SCIP_CALL( SCIPsolve(subscip) );

   /* terminate if the user has interupted the pricing problem */
   if ( SCIPpressedCtrlC(subscip) )
      return SCIP_OKAY;

   /* terminate if no violated inequality has been found */
   if ( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD )
   {
      SCIPdebugMsg(scip, "inequality is not violated -> no variable added\n");

      return SCIP_OKAY;
   }

   /* extract solution from subscip instance */
   subscipsol = SCIPgetBestSol(subscip);
   assert( subscipsol != NULL );

   subscipobj = SCIPgetSolOrigObj(subscip, subscipsol);

   if ( SCIPisLE(subscip, subscipobj, 1.0) )
   {
      SCIPdebugMsg(scip, "inequality is not violated -> no variable added\n");

      return SCIP_OKAY;
   }

   /* terminate if no violated inequality has been found */
   if ( SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD
      || SCIPisLE(subscip, subscipobj, 1.0) )
   {
      SCIPdebugMsg(scip, "inequality is not violated -> no variable added\n");

      return SCIP_OKAY;
   }

   /* extract the separated set from the solution and store it in standard format */
   SCIP_CALL( SCIPallocBufferArray(subscip, &separatedineq, dimension + 1) );
   separatedineq[0] = SCIPgetSolVal(subscip, subscipsol, ineqvars[dimension]);
   for (i = 0; i < dimension; ++i)
      separatedineq[i+1] = -SCIPgetSolVal(subscip, subscipsol, ineqvars[i]);

   SCIP_CALL( SCIPallocBufferArray(scip, &consids, nY) );
   nconss = 0;

   /* check which points in Y violate this inequality */
   for (i = 0; i < nY; ++i)
   {
      val = separatedineq[0];
      for (j = 0; j < dimension; ++j)
         val += separatedineq[j+1] * Y->points[i][j];

      /* stop, inequality is not violated */
      if ( SCIPisGT(subscip, val, -eps) )
         continue;

      consids[nconss++] = i;
   }

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   /* possibly search for a maximally separable set */
   if ( maxsepa )
   {
      SCIP_Bool success;
      SCIP_Longint nodelimit;
      SCIP_Real* inequ;
      SCIP_Real timelimit;
      int* fixed;
      SCIP_Real* fixedvals;
      int nfixed = 0;
      int* viol;
      int nviol;

      SCIP_CALL( SCIPallocBufferArray(scip, &fixed, Y->ndatapoints) );
      SCIP_CALL( SCIPallocBufferArray(scip, &fixedvals, Y->ndatapoints) );

      for (i = 0; i < Y->ndatapoints; ++i)
      {
         val = separatedineq[0];
         for (j = 0; j < dimension; ++j)
            val += separatedineq[j+1] * Y->points[i][j];

         /* stop, inequality is not violated */
         if ( SCIPisGT(subscip, val, -eps) )
            continue;

         fixed[nfixed] = i;
         fixedvals[nfixed++] = 1.0;
      }

      SCIP_CALL( SCIPgetLongintParam(scip, "rc/pricerimprovenodelimit", &nodelimit) );
      SCIP_CALL( SCIPgetRealParam(scip, "limits/time", &timelimit) );
      timelimit -= SCIPgetSolvingTime(scip);

      if ( SCIPisGT(scip, timelimit, 0.0) )
      {
         SCIP_CALL( SCIPallocBufferArray(scip, &inequ, X->dimension + 1) );
         SCIP_CALL( SCIPallocBufferArray(scip, &viol, Y->ndatapoints) );
         SCIP_CALL( SCIPfindMaximalSeparatingInequality(scip, X, Y, X->ndatapoints, Y->ndatapoints,
               NULL, X->dimension, absmaxX, fixed, fixedvals, nfixed, viol, &nviol, inequ, conshdlr,
               nodelimit, timelimit, &success) );

         /* we should at least be able to reproduce the already found inequality */
         if ( success )
         {
            for (i = 0; i <= X->dimension; ++i)
               separatedineq[i] = inequ[i];
            for (i = 0; i < nviol; ++i)
               consids[i] = viol[i];
            nconss = nviol;
         }

         SCIPfreeBufferArray(scip, &viol);
         SCIPfreeBufferArray(scip, &inequ);
      }

      SCIPfreeBufferArray(scip, &fixedvals);
      SCIPfreeBufferArray(scip, &fixed);
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "patternvar_%d", SCIPprobdataGetNVars(probdata));

   SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, consids, nconss, separatedineq, dimension + 1) );
   SCIP_CALL( SCIPcreateVarBinpacking(scip, &newvar, name, 1.0, FALSE, TRUE, TRUE, vardata) );

   SCIP_CALL( SCIPaddPricedVar(scip, newvar, 1.0) );
   SCIP_CALL( SCIPchgVarUbLazy(scip, newvar, 1.0) );

   *addvar = TRUE;

   /* adapt master problem */
   for (i = 0; i < nconss; ++i)
   {
      if( ! SCIPconsIsEnabled(coverconss[i]) )
         continue;

      SCIP_CALL( SCIPaddCoefSetppc(scip, coverconss[consids[i]], newvar) );
   }

   SCIP_CALL( SCIPreleaseVar(scip, &newvar) );

   SCIPfreeBufferArray(scip, &consids);
   SCIPfreeBufferArray(subscip, &separatedineq);

   return SCIP_OKAY;
}


/*
 * Callback methods of variable pricer
 */

/** destructor of variable pricer to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PRICERFREE(pricerFreePattern)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;

   assert( scip != NULL );
   assert( pricer != NULL );

   pricerdata = SCIPpricerGetData(pricer);

   if ( pricerdata != NULL )
   {
      /* free memory */
      SCIPfreeBlockMemoryArrayNull(scip, &pricerdata->coverconss, pricerdata->ncoverconss);

      SCIPfreeBlockMemory(scip, &pricerdata);
   }

   return SCIP_OKAY;
}


/** initialization method of variable pricer (called after problem was transformed) */
static
SCIP_DECL_PRICERINIT(pricerInitPattern)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   SCIP_CONS* cons;
   int c;

   assert( scip != NULL );
   assert( pricer != NULL );

   pricerdata = SCIPpricerGetData(pricer);
   assert( pricerdata != NULL );

   /* get transformed constraints */
   for (c = 0; c < pricerdata->ncoverconss; ++c)
   {
      cons = pricerdata->coverconss[c];

      /* release original constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &pricerdata->coverconss[c]) );

      /* get transformed constraint */
      SCIP_CALL( SCIPgetTransformedCons(scip, cons, &pricerdata->coverconss[c]) );

      /* capture transformed constraint */
      SCIP_CALL( SCIPcaptureCons(scip, pricerdata->coverconss[c]) );
   }

   return SCIP_OKAY;
}


/** solving process deinitialization method of variable pricer (called before branch and bound process data is freed) */
static
SCIP_DECL_PRICEREXITSOL(pricerExitsolPattern)
{  /*lint --e{715}*/
   SCIP_PRICERDATA* pricerdata;
   int c;

   assert( scip != NULL );
   assert( pricer != NULL );

   pricerdata = SCIPpricerGetData(pricer);
   assert( pricerdata != NULL );

   /* get release constraints */
   for (c = 0; c < pricerdata->ncoverconss; ++c)
   {
      /* release constraint */
      SCIP_CALL( SCIPreleaseCons(scip, &(pricerdata->coverconss[c])) );
   }

   return SCIP_OKAY;
}


/** reduced cost pricing method of variable pricer for feasible LPs */
static
SCIP_DECL_PRICERREDCOST(pricerRedcostPattern)
{  /*lint --e{715}*/
   SCIP* subscip;
   SCIP_PRICERDATA* pricerdata;
   SCIP_Bool addvar;
   SCIP_Bool atrootnode;
   SCIP_Bool success;

   assert( scip != NULL );
   assert( pricer != NULL );

   SCIPdebugMsg(scip, "start pricing round of pricer %s\n", PRICER_NAME);

   (*result) = SCIP_DIDNOTRUN;

   pricerdata = SCIPpricerGetData(pricer);
   assert( pricerdata != NULL );

   SCIP_CALL( SCIPcreate(&pricerdata->subscip) );
   subscip = pricerdata->subscip;

   atrootnode = SCIPnodeGetDepth(SCIPgetCurrentNode(scip)) == 0;

   SCIP_CALL( createPricingModel(scip, subscip, pricerdata, atrootnode) );
   SCIP_CALL( changeObjectivePricingProblem(scip, subscip, pricerdata->coverconss, pricerdata->nY, pricerdata->selectvars) );
   SCIP_CALL( addBranchingDecisionConss(scip, subscip, pricerdata->selectvars, pricerdata->conshdlr) );
   SCIP_CALL( addFixedVarsConss(scip, subscip, pricerdata->selectvars, pricerdata->coverconss, pricerdata->Y->ndatapoints) );
   SCIP_CALL( SCIPsolvePricingProblem(scip, subscip, pricerdata->X, pricerdata->nX, pricerdata->Y, pricerdata->nY,
         pricerdata->dimension, pricerdata->absmaxX, pricerdata->ineqvars, pricerdata->coverconss, pricerdata->conshdlr,
         &addvar) );

   success = SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL;
   success = success || SCIPgetStatus(subscip) == SCIP_STATUS_BESTSOLLIMIT;
   success = success || SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE;
   success = success || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD;

   if( addvar || success )
      (*result) = SCIP_SUCCESS;

   SCIP_CALL( freePricingModel(scip, subscip, pricerdata, atrootnode) );

   return SCIP_OKAY;
}

/*
 * variable pricer specific interface methods
 */

/** creates the pattern variable pricer and includes it in SCIP */
SCIP_RETCODE SCIPincludePricerPattern(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PRICERDATA* pricerdata;
   SCIP_PRICER* pricer;

   /* create stable set variable pricer data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &pricerdata) );

   pricerdata->coverconss = NULL;
   pricerdata->ncoverconss = 0;
   pricerdata->nhidingsetidx = 0;
   pricerdata->maxnhidingsetidx = 0;
   pricerdata->nhidingsetcuts = 0;

   /* include variable pricer */
   SCIP_CALL( SCIPincludePricerBasic(scip, &pricer, PRICER_NAME, PRICER_DESC, PRICER_PRIORITY, PRICER_DELAY,
         pricerRedcostPattern, NULL, pricerdata) );
   assert(pricer != NULL);

   /* set non fundamental callbacks via setter functions */
   SCIP_CALL( SCIPsetPricerFree(scip, pricer, pricerFreePattern) );
   SCIP_CALL( SCIPsetPricerInit(scip, pricer, pricerInitPattern) );
   SCIP_CALL( SCIPsetPricerExitsol(scip, pricer, pricerExitsolPattern) );

   return SCIP_OKAY;
}


/** adds problem specific data to pricer and activates pricer */
SCIP_RETCODE SCIPpricerPatternActivate(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   absmaxX,            /**< maximum absolute value of coordinate in X */
   SCIP_CONS**           coverconss,         /**< covering constraints */
   int                   ncoverconss,        /**< number of covering constraints */
   SCIP_CONSHDLR*        conshdlr            /**< pointer to same/diff conshdlr */
   )
{
   SCIP_PRICER* pricer;
   SCIP_PRICERDATA* pricerdata;
   int c;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( absmaxX > 0 );
   assert( coverconss != NULL );
   assert( ncoverconss > 0 );
   assert( conshdlr != NULL );

   pricer = SCIPfindPricer(scip, PRICER_NAME);
   assert( pricer != NULL );

   pricerdata = SCIPpricerGetData(pricer);
   assert( pricerdata != NULL );

   /* copy arrays */
   SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &pricerdata->coverconss, coverconss, ncoverconss) );
   pricerdata->ncoverconss = ncoverconss;

   /* capture all constraints */
   for (c = 0; c < ncoverconss; ++c)
   {
      SCIP_CALL( SCIPcaptureCons(scip, coverconss[c]) );
   }

   /* create pricing model */
   pricerdata->X = X;
   pricerdata->Y = Y;
   pricerdata->nX = X->ndatapoints;
   pricerdata->nY = Y->ndatapoints;
   pricerdata->dimension = X->dimension;
   pricerdata->absmaxX = absmaxX;
   pricerdata->nhidingsetcuts = 0;
   pricerdata->conshdlr = conshdlr;

   /* activate pricer */
   SCIP_CALL( SCIPactivatePricer(scip, pricer) );

   return SCIP_OKAY;
}
