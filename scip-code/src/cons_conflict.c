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

/**@file   cons_conflict.c
 * @brief  Constraint handler to handle conflicts
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"

#include "cons_conflict.h"
#include "datapoints.h"

/* fundamental constraint handler properties */
#define CONSHDLR_NAME          "conflict"
#define CONSHDLR_DESC          "constraint handler to handle conflicts"
#define CONSHDLR_SEPAPRIORITY    +40100 /**< priority of the constraint handler for separation */
#define CONSHDLR_ENFOPRIORITY        -1 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY       -1 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_EAGERFREQ          100 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */
#define CONSHDLR_SEPAFREQ            -1 /**< frequency for separating cuts; zero means to separate only in the root node */
#define CONSHDLR_PROPFREQ            -1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_DELAYSEPA        FALSE /**< should separation method be delayed, if other separators found cuts? */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_PROP_TIMING      SCIP_PROPTIMING_BEFORELP /**< propagation timing mask of the constraint handler */

#define DEFAULT_TRANSFERCUTS      FALSE /**< whether separated cuts for one class shall be transferred to other classes */

/*
 * Data structures
 */

/** constraint data for conflict constraints */
struct SCIP_ConsData
{
   Datapoints*           X;                  /**< feasible points */
   Datapoints*           Y;                  /**< infeasible points */
   int                   nY;                 /**< number of infeasible points */
   int                   nclasses;           /**< number of classes of points */
   int                   absmax;             /**< maximum absolute coordinate entry of a point in X */
   SCIP_VAR***           violvars;           /**< (ninfeasible x nclasses) matrix of violation variables */
};

/** constraint handler data */
struct SCIP_ConshdlrData
{
   SCIP_Bool             transfercuts;       /**< whether separated cuts for one class shall be transferred to other classes */
};


/*
 * Local methods
 */


/** frees an orbitope constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to orbitope constraint data */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}


/** creates conflict constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONSDATA**       consdata,           /**< pointer to constraint data */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   nclasses,           /**< number of classes of points */
   SCIP_VAR***           violvars,           /**< (ninfeasible x nclasses) matrix of violation variables */
   int                   absmax              /**< maximum absolute coordinate entry of a point in X */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nclasses > 0 );
   assert( violvars != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->X = X;
   (*consdata)->Y = Y;
   (*consdata)->nY = Y->ndatapoints;
   (*consdata)->nclasses = nclasses;
   (*consdata)->absmax = absmax;
   (*consdata)->violvars = violvars;

   return SCIP_OKAY;
}


/** performs the check of separability */
static
SCIP_RETCODE doIsSeparable(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int*                  selectedpoints,     /**< array of selected infeasible points */
   int                   nselectedpoints,    /**< number of selected infeasible points */
   int                   absmax,             /**< maximum absolute coordinate entry of a point in X */
   SCIP_Real             eps,                /**< epsilon used for checking separability */
   SCIP_Bool*            separable           /**< pointer to store whether solution is separable */
   )
{
   SCIP* subscip;
   SCIP_VAR** ineqvars;
   SCIP_CONS** validconss;
   SCIP_CONS** sepaconss;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nselectedpoints > 0 );
   assert( separable != NULL );

   SCIP_CALL( SCIPallocBufferArray(scip, &ineqvars, X->dimension + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &validconss, X->ndatapoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaconss, nselectedpoints) );
   SCIP_CALL( SCIPcreate(&subscip) );

   *separable = FALSE;

   SCIP_CALL( SCIPcreateSeparationLP(subscip, X, Y, selectedpoints, nselectedpoints, eps, absmax,
         ineqvars, validconss, sepaconss) );

   SCIP_CALL( SCIPsolve(subscip) );

   if ( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      *separable = TRUE;
   else if ( ! (SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD) )
   {
      SCIPerrorMessage("solution of separation LP cannot be handled\n");
      SCIP_CALL( SCIPinterruptSolve(scip) );
      SCIPfreeBufferArray(scip, &sepaconss);
      SCIPfreeBufferArray(scip, &validconss);
      SCIPfreeBufferArray(scip, &ineqvars);
      return SCIP_OKAY;
   }

   SCIP_CALL( SCIPfreeSeparationLP(subscip, X->ndatapoints, nselectedpoints, X->dimension,
         ineqvars, validconss, sepaconss) );
   SCIP_CALL( SCIPfree(&subscip) );

   SCIPfreeBufferArray(scip, &sepaconss);
   SCIPfreeBufferArray(scip, &validconss);
   SCIPfreeBufferArray(scip, &ineqvars);

   return SCIP_OKAY;
}


/** sparsifies conflict cut */
static
SCIP_RETCODE sparsifyCut(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int*                  selectedpoints,     /**< array of selected infeasible points */
   int*                  nselectedpoints,    /**< pointer to store number of selected infeasible points */
   int                   absmax,             /**< maximum absolute coordinate entry of a point in X */
   SCIP_Real             eps                 /**< epsilon used for checking separability */
   )
{
   int removedidx;
   int i;
   SCIP_Bool separable = TRUE;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( selectedpoints != NULL );
   assert( nselectedpoints != NULL );
   assert( *nselectedpoints >= 2 );

   if ( *nselectedpoints == 2 )
      return SCIP_OKAY;

   /* select greedily a minimal set of points that is not separable */
   for (i = 2; i < *nselectedpoints; ++i)
   {
      SCIP_CALL( doIsSeparable(scip, X, Y, selectedpoints, i, absmax, eps, &separable) );

      if ( ! separable )
      {
         *nselectedpoints = i;
         break;
      }
   }

   /* sparsify cut further by removing points that do not affect non-separability;
    * skip the last point, because we know that we need it for non-separability
    */
   for (i = *nselectedpoints - 2; i >= 0; --i)
   {
      removedidx = selectedpoints[i];
      selectedpoints[i] = selectedpoints[*nselectedpoints - 1];

      SCIP_CALL( doIsSeparable(scip, X, Y, selectedpoints, *nselectedpoints - 1, absmax, eps, &separable) );

      /* add removeidx again because we need it for non-separability or reduce number of points */
      if ( separable )
         selectedpoints[*nselectedpoints - 1] = removedidx;
      else
         *nselectedpoints -= 1;
   }

   return SCIP_OKAY;
}

/** checks whether a set of infeasible points can be separated simultaneously from X */
static
SCIP_RETCODE isSeparable(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS*            cons,               /**< conflict constraint */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   classidx,           /**< the class for which the cut is violated */
   int                   nclasses,           /**< number of inequality classes */
   int*                  selectedpoints,     /**< array of selected infeasible points */
   int                   nselectedpoints,    /**< number of selected infeasible points */
   SCIP_VAR***           violvars,           /**< matrix of variables (or NULL if not needed) */
   int                   absmax,             /**< maximum absolute coordinate entry of a point in X */
   SCIP_Real             eps,                /**< epsilon used for checking separability */
   SCIP_Bool             checkonly,          /**< if not, also a cut is produced */
   SCIP_Bool*            solinfeasible,      /**< pointer to store whether solution is infeasible */
   SCIP_Bool*            nodeinfeasible      /**< pointer to store whether local subproblem is infeasible */
   )
{
   SCIP* subscip;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_VAR** ineqvars;
   SCIP_CONS** validconss;
   SCIP_CONS** sepaconss;
   SCIP_ROW* row;
   int nselectedpointsorig;
   int i;
   int j;
   SCIP_Bool transfercuts;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nclasses > 0 );
   assert( selectedpoints != NULL );
   assert( nselectedpoints >= 0 );
   assert( solinfeasible != NULL );
   assert( nodeinfeasible != NULL );

   *solinfeasible = FALSE;
   *nodeinfeasible = FALSE;
   nselectedpointsorig = nselectedpoints;

   /* treat trivial cases */
   if ( nselectedpoints <= 1 )
      return SCIP_OKAY;

   SCIP_CALL( SCIPallocBufferArray(scip, &ineqvars, X->dimension + 1) );
   SCIP_CALL( SCIPallocBufferArray(scip, &validconss, X->ndatapoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &sepaconss, nselectedpoints) );
   SCIP_CALL( SCIPcreate(&subscip) );

   SCIP_CALL( SCIPcreateSeparationLP(subscip, X, Y, selectedpoints, nselectedpoints, eps, absmax,
         ineqvars, validconss, sepaconss) );

   SCIP_CALL( SCIPsolve(subscip) );

   if ( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      goto FREEDATASTRUCTURES;
   else if ( ! (SCIPgetStatus(subscip) == SCIP_STATUS_INFEASIBLE || SCIPgetStatus(subscip) == SCIP_STATUS_INFORUNBD) )
   {
      SCIPerrorMessage("solution of separation LP cannot be handled\n");
      SCIP_CALL( SCIPinterruptSolve(scip) );
      goto FREEDATASTRUCTURES;
   }

   /* the problem is bounded, hence: selected points cannot be separated from X */
   *solinfeasible = TRUE;

   if ( checkonly )
      goto FREEDATASTRUCTURES;

   /* add conflict cut for all classes */
   assert( violvars != NULL );

   SCIP_CALL( sparsifyCut(scip, X, Y, selectedpoints, &nselectedpoints, absmax, eps) );

   conshdlr = SCIPconsGetHdlr(cons);
   assert( conshdlr != NULL );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   transfercuts = conshdlrdata->transfercuts;

   for (j = 0; j < nclasses; ++j)
   {
      if ( !transfercuts && j != classidx )
         continue;

      SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, cons, "conflict", -SCIPinfinity(scip), nselectedpoints - 1,
            FALSE, FALSE, TRUE) );
      SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

      for (i = 0; i < nselectedpoints; ++i)
      {
         SCIP_CALL( SCIPaddVarToRow(scip, row, violvars[selectedpoints[i]][j], 1.0) );
      }
      SCIP_CALL( SCIPflushRowExtensions(scip, row) );

      SCIP_CALL( SCIPaddRow(scip, row, FALSE, nodeinfeasible) );
      SCIP_CALL( SCIPreleaseRow(scip, &row) );
   }

 FREEDATASTRUCTURES:
   SCIP_CALL( SCIPfreeSeparationLP(subscip, X->ndatapoints, nselectedpointsorig, X->dimension,
         ineqvars, validconss, sepaconss) );
   SCIP_CALL( SCIPfree(&subscip) );

   SCIPfreeBufferArray(scip, &sepaconss);
   SCIPfreeBufferArray(scip, &validconss);
   SCIPfreeBufferArray(scip, &ineqvars);

   return SCIP_OKAY;
}


/**< check or enforce a given solution */
static
SCIP_RETCODE checkOrEnfo(
   SCIP*                 scip,               /**< SCIP instance */
   SCIP_CONS**           conss,              /**< list of conflict constraints */
   int                   nconss,             /**< number of conflict constraints */
   SCIP_SOL*             sol,                /**< solution to be checked/enforced */
   SCIP_Bool             check,              /**< whether we just check the solution */
   int*                  ngen,               /**< number of generated conflicts */
   SCIP_Bool*            infeasible          /**< whether the local node is infeasible */
   )
{
   SCIP_Real eps;
   int c;

   assert( scip != NULL );
   assert( conss != NULL );

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   for (c = 0; c < nconss; ++c)
   {
      SCIP_CONSDATA* consdata;
      SCIP_VAR*** violvars;
      int* selectedpoints;
      int nselectedpoints;
      int nY;
      int nclasses;
      int i;
      int j;
      SCIP_Bool solinfeasible = FALSE;
      SCIP_Bool nodeinfeasible = FALSE;

      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );
      assert( consdata->Y != NULL );
      assert( consdata->violvars != NULL );

      nY = consdata->nY;
      nclasses = consdata->nclasses;
      violvars = consdata->violvars;

      /* for each class, collect selected infeasible points */
      SCIP_CALL( SCIPallocBufferArray(scip, &selectedpoints, nY) );
      for (j = 0; j < nclasses && ! nodeinfeasible; ++j)
      {
         nselectedpoints = 0;

         for (i = 0; i < nY; ++i)
         {
            if ( SCIPgetSolVal(scip, sol, violvars[i][j]) > 0.5 )
               selectedpoints[nselectedpoints++] = i;
         }

         SCIP_CALL( isSeparable(scip, conss[c], consdata->X, consdata->Y, j, nclasses,
               selectedpoints, nselectedpoints, violvars, consdata->absmax, eps,
               check, &solinfeasible, &nodeinfeasible) );
         if ( nodeinfeasible )
            *infeasible = TRUE;
         else if ( solinfeasible )
            *ngen += 1;
      }

      SCIPfreeBufferArray(scip, &selectedpoints);

      if ( nodeinfeasible )
      {
         *infeasible = TRUE;
         break;
      }
   }

   return SCIP_OKAY;
}

/** propagation method for a single conflict constraint */
static
SCIP_RETCODE propagateCons(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons,               /**< constraint to be processed */
   SCIP_Bool*            infeasible,         /**< pointer to store TRUE, if the node can be cut off */
   int*                  nfixedvars          /**< pointer to add up the number of found domain reductions */
   )
{
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   Datapoints* X;
   Datapoints* Y;
   SCIP_Real eps;
   int* selectedpoints;
   int nselectedpoints;
   int nclasses;
   int nY;
   int c;
   int i;

   assert( scip != NULL );
   assert( cons != NULL );
   assert( infeasible != NULL );
   assert( nfixedvars != NULL );

   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );

   X = consdata->X;
   Y = consdata->Y;
   nY = consdata->nY;
   nclasses = consdata->nclasses;
   vars = consdata->violvars;
   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   *infeasible = FALSE;

   SCIP_CALL( SCIPallocBufferArray(scip, &selectedpoints, nY) );

   /* per class of inequalities, check whether the selected points cannot be separated already
    * or whether adding another point will lead to an inseparable set
    */
   for (c = 0; c < nclasses; ++c)
   {
      SCIP_Bool nodeinfeasible = FALSE;
      nselectedpoints = 0;

      for (i = 0; i < nY; ++i)
      {
         if ( SCIPvarGetLbLocal(vars[i][c]) > 0.5 )
            selectedpoints[nselectedpoints++] = i;
      }

      /* avoid trivial cases */
      if ( nselectedpoints <= 1 )
         continue;

      SCIP_CALL( isSeparable(scip, cons, X, Y, c, nclasses, selectedpoints, nselectedpoints, vars,
            consdata->absmax, eps, TRUE, infeasible, &nodeinfeasible) );

      /* stop, the 1-branched points cannot be separated */
      if ( *infeasible )
         break;

      /* check whether additional points can be fixed */
      for (i = 0; i < nY; ++i)
      {
         SCIP_Bool tmpinfeasible = FALSE;

         if ( SCIPvarGetLbLocal(vars[i][c]) < 0.5 && SCIPvarGetUbLocal(vars[i][c]) > 0.5 )
            selectedpoints[nselectedpoints] = i;
         else
            continue;

         SCIP_CALL( isSeparable(scip, cons, X, Y, c, nclasses, selectedpoints, nselectedpoints + 1, vars,
               consdata->absmax, eps, TRUE, &tmpinfeasible, &nodeinfeasible) );

         if ( tmpinfeasible )
         {
            SCIP_CALL( SCIPchgVarUb(scip, vars[i][c], 0.0) );
            *nfixedvars += 1;
         }
      }
   }

   SCIPfreeBufferArray(scip, &selectedpoints);

   return SCIP_OKAY;
}

/** separate constraints heuristically */
static
SCIP_RETCODE separateConstraints(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSHDLR*        conshdlr,           /**< constraint handler */
   SCIP_CONS**           conss,              /**< constraints to process */
   int                   nconss,             /**< number of constraints */
   SCIP_SOL*             sol,                /**< solution to separate (NULL for the LP solution) */
   SCIP_RESULT*          result              /**< pointer to store the result (should be initialized) */
   )
{
   SCIP_Bool infeasible = FALSE;
   int ncuts = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   /* loop through constraints */
   for (c = 0; c < nconss && ! infeasible; c++)
   {
      SCIP_CONSDATA* consdata;
      Datapoints* X;
      Datapoints* Y;
      SCIP_VAR*** vars;
      SCIP_Real* solution;
      int* pointorder;
      int* selectedpoints;
      SCIP_Real eps;
      int absmax;
      int nselectedpoints;
      int nY;
      int nclasses;
      int i;
      int j;
      SCIP_Real lhs;
      SCIP_Bool separable = TRUE;

      assert( conss[c] != NULL );

      /* get data of constraint */
      consdata = SCIPconsGetData(conss[c]);
      assert( consdata != NULL );

      X = consdata->X;
      Y = consdata->Y;
      nY = consdata->nY;
      nclasses = consdata->nclasses;
      vars = consdata->violvars;
      absmax = consdata->absmax;
      SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

      SCIP_CALL( SCIPallocBufferArray(scip, &pointorder, nY) );
      for (i = 0; i < nY; ++i)
         pointorder[i] = i;

      SCIP_CALL( SCIPallocBufferArray(scip, &selectedpoints, nY) );
      SCIP_CALL( SCIPallocBufferArray(scip, &solution, nY) );

      /* for each class, try to find a separating inequality */
      for (j = 0; j < nclasses && separable; ++j)
      {
         /* sort points in Y w.r.t. their value in the solution
          * use negative solution value to sort non-increasingly
          **/
         for (i = 0; i < nY; ++i)
            solution[i] = - SCIPgetSolVal(scip, sol, vars[pointorder[i]][j]);

         SCIPsortRealInt(solution, pointorder, nY);

         /* heuristically create a set of points that cannot be separated */
         lhs = 0.0;
         nselectedpoints = 0;
         for (i = 0; i < nY && separable; ++i)
         {
            lhs -= solution[i];

            /* terminate early if no violated inequality can exist */
            if ( SCIPisLE(scip, lhs, (SCIP_Real) i - 1) )
               break;

            /* add pointorder[i] to set of selected points and check for separability */
            selectedpoints[nselectedpoints++] = pointorder[i];
            SCIP_CALL( doIsSeparable(scip, X, Y, selectedpoints, nselectedpoints,
                  absmax, eps, &separable) );
         }

         /* possibly generate a cut */
         if ( ! separable )
         {
            SCIP_ROW* row;

            SCIP_CALL( SCIPcreateEmptyRowCons(scip, &row, conss[c], "conflictcut",
                  -SCIPinfinity(scip), (SCIP_Real) nselectedpoints - 1, FALSE, FALSE, TRUE) );
            SCIP_CALL( SCIPcacheRowExtensions(scip, row) );

            for (i = 0; i < nselectedpoints; ++i)
            {
               SCIP_CALL( SCIPaddVarToRow(scip, row, vars[pointorder[i]][j], 1.0) );
            }
            SCIP_CALL( SCIPflushRowExtensions(scip, row) );

            SCIP_CALL( SCIPaddRow(scip, row, FALSE, &infeasible) );
            SCIP_CALL( SCIPreleaseRow(scip, &row) );
            ++ncuts;
         }
      }

      SCIPfreeBufferArray(scip, &solution);
      SCIPfreeBufferArray(scip, &selectedpoints);
      SCIPfreeBufferArray(scip, &pointorder);
   }

   if ( infeasible )
   {
      SCIPdebugMsg(scip, "Infeasible node.\n");
      *result = SCIP_CUTOFF;
   }
   else if ( ncuts > 0 )
   {
      SCIPdebugMsg(scip, "Separated %dinequalities.\n", ncuts);
      *result = SCIP_SEPARATED;
   }
   else
   {
      SCIPdebugMsg(scip, "No violated inequality found during separation.\n");
   }

   return SCIP_OKAY;
}


/*
 * Callback methods of constraint handler
 */


/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteConflict)
{  /*lint --e{715}*/
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** frees constraint handler */
static
SCIP_DECL_CONSFREE(consFreeConflict)
{  /*lint --e{715}*/
   SCIP_CONSHDLRDATA* conshdlrdata;

   assert( scip != 0 );
   assert( conshdlr != 0 );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );

   conshdlrdata = SCIPconshdlrGetData(conshdlr);
   assert( conshdlrdata != NULL );

   SCIPfreeBlockMemory(scip, &conshdlrdata);

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransConflict)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create linear constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata, sourcedata->X, sourcedata->Y, sourcedata->nclasses,
         sourcedata->violvars, sourcedata->absmax) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for LP solutions */
static
SCIP_DECL_CONSENFOLP(consEnfolpConflict)
{  /*lint --e{715}*/
   int ngen = 0;
   SCIP_Bool infeasible = FALSE;

   *result = SCIP_FEASIBLE;

   SCIP_CALL( checkOrEnfo(scip, conss, nconss, NULL, FALSE, &ngen, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( ngen > 0 )
      *result = SCIP_SEPARATED;

   return SCIP_OKAY;
}


/** constraint enforcing method of constraint handler for pseudo solutions */
static
SCIP_DECL_CONSENFOPS(consEnfopsConflict)
{  /*lint --e{715}*/
   int ngen = 0;
   SCIP_Bool infeasible;
   *result = SCIP_FEASIBLE;

   SCIP_CALL( checkOrEnfo(scip, conss, nconss, NULL, TRUE, &ngen, &infeasible) );

   if ( ngen > 0 )
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** feasibility check method of constraint handler for integral solutions */
static
SCIP_DECL_CONSCHECK(consCheckConflict)
{  /*lint --e{715}*/
   int ngen = 0;
   SCIP_Bool infeasible;
   *result = SCIP_FEASIBLE;

   SCIP_CALL( checkOrEnfo(scip, conss, nconss, sol, TRUE, &ngen, &infeasible) );

   if ( ngen > 0 )
      *result = SCIP_INFEASIBLE;

   return SCIP_OKAY;
}


/** variable rounding lock method of constraint handler */
static
SCIP_DECL_CONSLOCK(consLockConflict)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_VAR*** vars;
   int nY;
   int nclasses;
   int i;
   int j;

   assert( scip != NULL );

   consdata = SCIPconsGetData(cons);
   assert( consdata != NULL );
   assert( consdata->Y != NULL );

   nY = consdata->nY;
   nclasses = consdata->nclasses;
   vars = consdata->violvars;

   assert( nY > 0 );
   assert( nclasses > 0 );
   assert( vars != NULL );

   for (i = 0; i < nY; ++i)
   {
      for (j = 0; j < nclasses; ++j)
      {
         SCIP_CALL( SCIPaddVarLocksType(scip, vars[i][j], locktype, nlocksneg, nlockspos) );
      }
   }

   return SCIP_OKAY;
}


/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropConflict)
{  /*lint --e{715}*/
   SCIP_Bool infeasible = FALSE;
   int nfixedvars = 0;
   int c;

   assert( scip != NULL );
   assert( conshdlr != NULL );
   assert( strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* propagate all constraints */
   for (c = 0; c < nconss && !infeasible; ++c)
   {
      assert( conss[c] != 0 );

      SCIPdebugMsg(scip, "Propagation of conflict constraint <%s> ...\n", SCIPconsGetName(conss[c]));

      SCIP_CALL( propagateCons(scip, conss[c], &infeasible, &nfixedvars) );
   }

   /* return the correct result */
   if ( infeasible )
   {
      *result = SCIP_CUTOFF;
      SCIPdebugMsg(scip, "Propagation proved node to be infeasible.\n");
   }
   else if ( nfixedvars > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      SCIPdebugMsg(scip, "Propagated %d variables.\n", nfixedvars);
   }
   else if ( nconss > 0 )
   {
      *result = SCIP_DIDNOTFIND;
      SCIPdebugMsg(scip, "Propagation did not find anything.\n");
   }

   return SCIP_OKAY;
}


/** separation method of constraint handler for LP solutions */
static
SCIP_DECL_CONSSEPALP(consSepalpConflict)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of orbitope constraint handler <%s> for LP solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTRUN;

   /* if solution is integer, skip separation */
   if ( SCIPgetNLPBranchCands(scip) <= 0 )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, NULL, result) );

   return SCIP_OKAY;
}

/** separation method of constraint handler for arbitrary primal solutions */
static
SCIP_DECL_CONSSEPASOL(consSepasolConflict)
{  /*lint --e{715}*/
   assert( scip != NULL );
   assert( result != NULL );

   SCIPdebugMsg(scip, "Separation of conflict constraint handler <%s> for primal solution.\n", SCIPconshdlrGetName(conshdlr));

   *result = SCIP_DIDNOTFIND;

   /* separate constraints */
   SCIP_CALL( separateConstraints(scip, conshdlr, conss, nconss, sol, result) );

   return SCIP_OKAY;
}


/*
 * constraint specific interface methods
 */

/** creates the handler for conflict constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrConflict(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata;
   SCIP_CONSHDLR* conshdlr;

   /* create conflict constraint handler data */
   SCIP_CALL( SCIPallocBlockMemory(scip, &conshdlrdata) );

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpConflict, consEnfopsConflict, consCheckConflict, consLockConflict,
         conshdlrdata) );
   assert(conshdlr != NULL);

   /* set non-fundamental callbacks via specific setter functions */
   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteConflict) );
   SCIP_CALL( SCIPsetConshdlrFree(scip, conshdlr, consFreeConflict) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransConflict) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropConflict, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrSepa(scip, conshdlr, consSepalpConflict, consSepasolConflict, CONSHDLR_SEPAFREQ,
         CONSHDLR_SEPAPRIORITY, CONSHDLR_DELAYSEPA) );

   SCIP_CALL( SCIPaddBoolParam(scip, "constraints/conflict/transfercuts",
         "whether separated cuts for one class shall be transferred to other classes",
         &conshdlrdata->transfercuts, TRUE, DEFAULT_TRANSFERCUTS, NULL, NULL) );

   return SCIP_OKAY;
}

/** creates and captures a conflict constraint
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   nclasses,           /**< number of classes of points */
   SCIP_VAR***           violvars,           /**< (ninfeasible x nclasses) matrix of violation variables */
   int                   absmax,             /**< maximum absolute coordinate entry of a point in X */
   SCIP_Bool             initial,            /**< should the LP relaxation of constraint be in the initial LP?
                                              *   Usually set to TRUE. Set to FALSE for 'lazy constraints'. */
   SCIP_Bool             separate,           /**< should the constraint be separated during LP processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             enforce,            /**< should the constraint be enforced during node processing?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             check,              /**< should the constraint be checked for feasibility?
                                              *   TRUE for model constraints, FALSE for additional, redundant constraints. */
   SCIP_Bool             propagate,          /**< should the constraint be propagated during node processing?
                                              *   Usually set to TRUE. */
   SCIP_Bool             local,              /**< is constraint only valid locally?
                                              *   Usually set to FALSE. Has to be set to TRUE, e.g., for branching constraints. */
   SCIP_Bool             modifiable,         /**< is constraint modifiable (subject to column generation)?
                                              *   Usually set to FALSE. In column generation applications, set to TRUE if pricing
                                              *   adds coefficients to this constraint. */
   SCIP_Bool             dynamic,            /**< is constraint subject to aging?
                                              *   Usually set to FALSE. Set to TRUE for own cuts which
                                              *   are separated as constraints. */
   SCIP_Bool             removable,          /**< should the relaxation be removed from the LP due to aging or cleanup?
                                              *   Usually set to FALSE. Set to TRUE for 'lazy constraints' and 'user cuts'. */
   SCIP_Bool             stickingatnode      /**< should the constraint always be kept at the node where it was added, even
                                              *   if it may be moved to a more global node?
                                              *   Usually set to FALSE. Set to TRUE to for constraints that represent node data. */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the conflict constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if ( conshdlr == NULL )
   {
      SCIPerrorMessage("conflict constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create constraint data */
   consdata = NULL;
   SCIP_CALL( consdataCreate(scip, &consdata, X, Y, nclasses, violvars, absmax) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, initial, separate, enforce, check, propagate,
         local, modifiable, dynamic, removable, stickingatnode) );

   return SCIP_OKAY;
}

/** creates and captures a conflict constraint with all its constraint flags set to their
 *  default values
 *
 *  @note the constraint gets captured, hence at one point you have to release it using the method SCIPreleaseCons()
 */
SCIP_RETCODE SCIPcreateConsBasicConflict(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int                   nclasses,           /**< number of classes of points */
   SCIP_VAR***           violvars,           /**< (ninfeasible x nclasses) matrix of violation variables */
   int                   absmax              /**< maximum absolute coordinate entry of a point in X */
   )
{
   SCIP_CALL( SCIPcreateConsConflict(scip, cons, name, X, Y, nclasses, violvars, absmax,
         TRUE, TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE) );

   return SCIP_OKAY;
}

/** creates LP to check whether two sets of points can be separated */
SCIP_RETCODE SCIPcreateSeparationLP(
   SCIP*                 subscip,            /**< subSCIP instance to be created */
   Datapoints*           X,                  /**< feasible points */
   Datapoints*           Y,                  /**< infeasible points */
   int*                  selectedpoints,     /**< selected infeasible points */
   int                   nselectedpoints,    /**< number of selected infeasible points */
   SCIP_Real             eps,                /**< epsilon value used to classify separated points */
   int                   absmax,             /**< maximum absolute value of a coordinate in X */
   SCIP_VAR**            ineqvars,           /**< allocated array to store inequality variables */
   SCIP_CONS**           validconss,         /**< allocated array to store validity constraints */
   SCIP_CONS**           sepaconss           /**< allocated array to store separability constraints */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP_Real* vals;
   int dim;
   int nX;
   int i;
   int j;

   assert( X != NULL );
   assert( Y != NULL );
   assert( ineqvars != NULL );
   assert( validconss != NULL );
   assert( sepaconss != NULL );

   /* create subSCIP */
   SCIP_CALL( SCIPcreateProbBasic(subscip, "checkSeparability") );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* load plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   dim = X->dimension;
   nX = X->ndatapoints;

   assert( dim > 0 );
   assert( nX > 0 );

   /* create variables */
   for (i = 0; i < dim; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "lhs_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &ineqvars[i], name, -1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, ineqvars[i]) );
   }
   SCIP_CALL( SCIPcreateVar(subscip, &ineqvars[dim], "rhs", -dim * absmax, dim * absmax, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, ineqvars[dim]) );

   /* create constraints */
   SCIP_CALL( SCIPallocBufferArray(subscip, &vals, dim + 1) );
   vals[dim] = -1.0;

   /* create validconss */
   for (i = 0; i < nX; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "validcons_%d", i);

      for (j = 0; j < dim; ++j)
         vals[j] = X->points[i][j];

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &validconss[i], name,
            dim + 1, ineqvars, vals, -SCIPinfinity(subscip), 0.0) );
      SCIP_CALL( SCIPaddCons(subscip, validconss[i]) );
   }

   /* create sepaconss */
   vals[dim] = 1.0;
   for (i = 0; i < nselectedpoints; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sepacons_%d", i);

      for (j = 0; j < dim; ++j)
         vals[j] = -Y->points[selectedpoints[i]][j];

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &sepaconss[i], name,
            dim + 1, ineqvars, vals, -SCIPinfinity(subscip), -eps) );
      SCIP_CALL( SCIPaddCons(subscip, sepaconss[i]) );
   }

   SCIPfreeBufferArray(subscip, &vals);

   return SCIP_OKAY;
}


/** frees separation LP */
SCIP_RETCODE SCIPfreeSeparationLP(
   SCIP*                 subscip,            /**< subSCIP instance to be freed */
   int                   nfeaspoints,        /**< number of feasible points */
   int                   ninfpoints,         /**< number of selected infeasible points */
   int                   dim,                /**< dimension of problem */
   SCIP_VAR**            ineqvars,           /**< array of variables */
   SCIP_CONS**           validconss,         /**< array of valid conss */
   SCIP_CONS**           sepaconss           /**< array of separation conss */
   )
{
   int i;

   assert( subscip != NULL );
   assert( ineqvars != NULL );
   assert( validconss != NULL );
   assert( sepaconss != NULL );

   for (i = 0; i < nfeaspoints; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &validconss[i]) );
   }
   for (i = 0; i < ninfpoints; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(subscip, &sepaconss[i]) );
   }
   for (i = 0; i <= dim; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &ineqvars[i]) );
   }
   SCIP_CALL( SCIPfreeTransform(subscip) );

   return SCIP_OKAY;
}
