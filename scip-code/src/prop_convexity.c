/**@file   prop_convexity.c
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"
#include "math.h"
#include "convex_hull.h"

#include "datapoints.h"
#include "probdata_rc_compact.h"
#include "probdata_rc_conflict.h"
#include "prop_convexity.h"
#include "vardata_compact.h"
#include "typedefs.h"

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

/* fundamental propagator properties */
#define PROP_NAME              "convexity"
#define PROP_DESC              "propagates that sets of separated infeasible points are 'convex'"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                     1 /**< propagator frequency */
#define PROP_DELAY                FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define PROP_TIMING             SCIP_PROPTIMING_BEFORELP/**< propagation timing mask */

/* optional propagator properties */
#define PROP_PRESOL_PRIORITY         -1 /**< priority of the presolving method (>= 0: before, < 0: after constraint handlers); combined with presolvers */
#define PROP_PRESOLTIMING       SCIP_PRESOLTIMING_MEDIUM /* timing of the presolving method (fast, medium, or exhaustive) */
#define PROP_PRESOL_MAXROUNDS        -1 /**< maximal number of presolving rounds the presolver participates in (-1: no
                                         *   limit) */

/*
 * Data structures
 */

/** propagator data */
struct SCIP_PropData
{
   Datapoints*           X;                  /**< pointer to feasible points */
   Datapoints*           Y;                  /**< pointer to infeasible points */
   int                   nX;                 /**< number of feasible points */
   int                   nY;                 /**< number of infeasible points */
   int                   dimension;          /**< dimension of points */
   int                   ninequalities;      /**< number of inequalities in formualtion */
   SCIP_VAR***           vars;               /**< (nY x ninequalities)-array of variables */
   SCIP_NODE*            lastnode;           /**< last node considered in propagation */
   SCIP_Bool             enabled;            /**< whether the propagator is enabled */
};


/*
 * Local methods
 */

/** frees data of propagator */
static
SCIP_RETCODE propdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to propdata */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   SCIPfreeBlockMemory(scip, propdata);

   return SCIP_OKAY;
}


/** creates data structure of propagator */
static
SCIP_RETCODE propdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA**       propdata            /**< pointer to store propagator data */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );

   SCIP_CALL( SCIPallocBlockMemory(scip, propdata) );

   return SCIP_OKAY;
}


/** set data structure of propagator */
static
SCIP_RETCODE propdataSet(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< pointer to store propagator data */
   Datapoints*           X,                  /**< pointer to feasible points */
   Datapoints*           Y,                  /**< pointer to infeasible points */
   int                   nX,                 /**< number of feasible points */
   int                   nY,                 /**< number of infeasible points */
   int                   dimension,          /**< dimension of points */
   int                   ninequalities,      /**< number of inequalities in formualtion */
   SCIP_VAR***           vars,               /**< (nY x ninequalities)-array of variables */
   SCIP_NODE*            lastnode            /**< last node considered in propagation */
   )
{
   assert( scip != NULL );
   assert( propdata != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( dimension > 0 );
   assert( vars != NULL );

   propdata->X = X;
   propdata->Y = Y;
   propdata->nX = nX;
   propdata->nY = nY;
   propdata->dimension = dimension;
   propdata->ninequalities = ninequalities;
   propdata->vars = vars;
   propdata->lastnode = lastnode;

   return SCIP_OKAY;
}


/** performs the propagation based on intersection of convex hull of infeasible points */
static
SCIP_RETCODE propagate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROPDATA*        propdata,           /**< data of propagator */
   int*                  nfixings,           /**< pointer to store number of fixings found by propagator */
   SCIP_Bool*            infeasible          /**< pointer to store if infeasibility has been detected */
   )
{
   Datapoints* Y;
   int nY;
   SCIP_VAR*** vars;

   int* separatedpoints;
   int nseparatedpoints;
   int idx;
   SCIP_VAR* var;

   SCIP_NODE* node;
   SCIP_BOUNDCHG* boundchg;
   SCIP_DOMCHG* domchg;
   SCIP_VARDATA* vardata;
   SCIP_VAR* branchvar;
   int nboundchgs;
   int i;
   int f;
   int j;

   dd_MatrixPtr generators;
   dd_MatrixPtr facets;
   SCIP_Real val;
   SCIP_Bool success = TRUE;

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nfixings != NULL );
   assert( infeasible != NULL );

   *nfixings = 0;
   *infeasible = FALSE;

   Y = propdata->Y;
   nY = propdata->nY;
   vars = propdata->vars;

   assert( Y != NULL );
   assert( nY > 0 );
   assert( vars != NULL );

   /* do nothing if we are in a probing node */
   if ( SCIPinProbing(scip) )
      return SCIP_OKAY;

   /* get last branching decision */
   node = SCIPgetCurrentNode(scip);
   assert( node != NULL );

   domchg = SCIPnodeGetDomchg(node);
   if ( domchg == NULL )
      return SCIP_OKAY;

   nboundchgs = SCIPdomchgGetNBoundchgs(domchg);
   if ( nboundchgs == 0 )
      return SCIP_OKAY;

   boundchg = SCIPdomchgGetBoundchg(domchg, 0);

   /* branching decisions have to be in the beginning of the bound change array */
   if ( SCIPboundchgGetBoundchgtype(boundchg) != SCIP_BOUNDCHGTYPE_BRANCHING )
      return SCIP_OKAY;

   branchvar = SCIPboundchgGetVar(boundchg);

   /* skip non-binary branching variables */
   if ( SCIPvarGetType(branchvar) != SCIP_VARTYPE_BINARY )
      return SCIP_OKAY;

   vardata = SCIPvarGetData(branchvar);
   if ( vardata == NULL )
      return SCIP_OKAY;

   /* get inequality index of branched variable */
   idx = SCIPvardataGetInequalityidx(vardata);

   /* iterate over infeasible points and collect the ones to be separated by inequality idx */
   nseparatedpoints = 0;
   SCIP_CALL( SCIPallocBufferArray(scip, &separatedpoints, nY) );

   for (i = 0; i < nY; ++i)
   {
      var = vars[i][idx];

      if ( SCIPvarGetLbLocal(var) > 0.5 )
         separatedpoints[nseparatedpoints++] = i;
   }

   /* we cannot propagate if the set of separated points in not full-dimensional */
   if ( nseparatedpoints <= Y->dimension )
   {
      SCIPfreeBufferArray(scip, &separatedpoints);
      return SCIP_OKAY;
   }

   /* compute convex hull of the separated points */
   dd_set_global_constants();

   generators = constructGeneratorMatrixPoints(Y, separatedpoints, nseparatedpoints);
   facets = computeConvexHullFacets(scip, generators, &success);

   dd_FreeMatrix(generators);
   SCIPfreeBufferArray(scip, &separatedpoints);

   /* do not propagate if we could not compute an outer description of the convex hull,
    * or the convex hull is not full-dimensional
    */
   if ( ! success )
   {
      dd_free_global_constants();

      return SCIP_OKAY;
   }

   /* iterate over infeasible points not (yet) separated by inequality idx;
    * if they are contained in the convex hull, they are also separated by the inequality
    */
   for (i = 0; i < nY; ++i)
   {
      var = vars[i][idx];

      /* skip points already separated by idx */
      if ( SCIPvarGetLbLocal(var) > 0.5 )
         continue;

      /* check whether point is contained in convex hull */
      for (f = 0; f < facets->rowsize; ++f)
      {
         val = getReal(facets->matrix[f][0]);
         for (j = 1; j < facets->colsize; ++j)
            val += getReal(facets->matrix[f][j]) * Y->points[i][j-1];

         if ( SCIPisSumLT(scip, val, 0.0) )
            break;
      }

      /* skip points not contained in the convex hull */
      if ( f < facets->rowsize )
         continue;

      /* assign point to be separated by idx */
      if ( SCIPvarGetUbLocal(var) < 0.5 )
      {
         *infeasible = TRUE;
         break;
      }
      else if ( SCIPvarGetUbLocal(var) > 0.5 && SCIPvarGetLbLocal(var) < 0.5 )
      {
         SCIP_CALL( SCIPchgVarLb(scip, var, 1.0) );
         *nfixings += 1;
      }
   }

   dd_FreeMatrix(facets);
   dd_free_global_constants();

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitConvexity)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_PROPDATA* propdata;
   Datapoints* X;
   Datapoints* Y;
   int dimension;
   int ninequalities;
   SCIP_VAR*** vars;
   int method;

   probdata = SCIPgetProbData(scip);
   assert( probdata != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( SCIPgetBoolParam(scip, "propagating/convexity/enabled", &propdata->enabled) );

   if ( ! propdata->enabled )
      return SCIP_OKAY;

   SCIP_CALL( SCIPgetIntParam(scip, "rc/method", &method) );
   if ( method == METHOD_CONFLICT || method == METHOD_HYBRID_CONFLICT )
   {
      X = SCIPprobdataConflictGetX(probdata);
      Y = SCIPprobdataConflictGetY(probdata);
      dimension = X->dimension;
      ninequalities = SCIPprobdataConflictGetUb(probdata);
      vars = SCIPprobdataConflictGetViolatedvars(probdata);
   }
   else
   {
      X = SCIPprobdataCompactGetX(probdata);
      Y = SCIPprobdataCompactGetY(probdata);
      dimension = X->dimension;
      ninequalities = SCIPprobdataGetUb(probdata);
      vars = SCIPprobdataGetViolatedvars(probdata);
   }

   /* create convexity propagator data */
   SCIP_CALL( propdataSet(scip, propdata, X, Y, X->ndatapoints, Y->ndatapoints,
         dimension, ninequalities, vars, NULL) );

   return SCIP_OKAY;
}


/** destructor of propagator to free user data (called when SCIP is exiting) */
static
SCIP_DECL_PROPFREE(propFreeConvexity)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;

   assert( scip != NULL );
   assert( prop != NULL );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   SCIP_CALL( propdataFree(scip, &propdata) );

   return SCIP_OKAY;
}


/** presolving method of propagator */
static
SCIP_DECL_PROPPRESOL(propPresolConvexity)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   int nfixings = 0;
   SCIP_Bool infeasible = FALSE;

   assert( scip != NULL );
   assert( prop != NULL );
   assert( result != NULL );
   assert( SCIPgetStage(scip) == SCIP_STAGE_PRESOLVING );

   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   *result = SCIP_DIDNOTRUN;

   if ( ! propdata->enabled )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nfixings > 0 )
   {
      *result = SCIP_REDUCEDDOM;
      *nfixedvars += nfixings;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecConvexity)
{  /*lint --e{715}*/
   SCIP_PROPDATA* propdata;
   SCIP_Bool infeasible = FALSE;
   int nfixings = 0;

   assert( scip != NULL );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   /* do not run if we are in the root or not yet solving */
   if ( SCIPgetDepth(scip) <= 0 || SCIPgetStage(scip) < SCIP_STAGE_SOLVING )
      return SCIP_OKAY;

   /* get data */
   propdata = SCIPpropGetData(prop);
   assert( propdata != NULL );

   if ( ! propdata->enabled )
      return SCIP_OKAY;

   *result = SCIP_DIDNOTFIND;

   SCIP_CALL( propagate(scip, propdata, &nfixings, &infeasible) );

   if ( infeasible )
      *result = SCIP_CUTOFF;
   else if ( nfixings > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** propagation conflict resolving method of propagator */
static
SCIP_DECL_PROPRESPROP(propRespropConvexity)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the intersection propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropConvexity(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( propdataCreate(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecConvexity, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitConvexity) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeConvexity) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolConvexity, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropConvexity) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/convexity/enabled",
         "whether the proagator is enabled", NULL, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
