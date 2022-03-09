/**@file   prop_intersection.c
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
#include "prop_intersection.h"
#include "typedefs.h"
#include "vardata_compact.h"

#include "scip/cons_linear.h"
#include "scip/scipdefplugins.h"

/* fundamental propagator properties */
#define PROP_NAME              "intersection"
#define PROP_DESC              "propagates information on infeasible points based on intersections"
#define PROP_PRIORITY                 0 /**< propagator priority */
#define PROP_FREQ                    -1 /**< propagator frequency */
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


/** checks whether two full-dimensional polyhedra in H-representation intersect */
static
SCIP_Bool polyhedraIntersect(
   dd_MatrixPtr          facets1,            /**< H-description of first polyhedron */
   dd_MatrixPtr          facets2             /**< H-description of second polyhedron */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP* subscip;
   dd_rowrange rowidx;
   dd_colrange colidx;
   SCIP_VAR** vars;
   SCIP_Real* coefs;
   SCIP_CONS* cons;
   dd_colrange nvars;
   SCIP_Bool intersect = FALSE;

   nvars = facets1->colsize - 1;

   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &vars, nvars) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &coefs, nvars) );

   SCIP_CALL( SCIPcreateProbBasic(subscip, "intersect") );

   /* add variables */
   for (colidx = 0; colidx < nvars; ++colidx)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "var_%d", colidx);

      SCIP_CALL( SCIPcreateVar(subscip, &vars[colidx], name, -SCIPinfinity(subscip), SCIPinfinity(subscip), 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, vars[colidx]) );
   }

   /* generate constraints */
   for (rowidx = 0; rowidx < facets1->rowsize; ++rowidx)
   {
      for (colidx = 0; colidx < nvars; ++colidx)
         coefs[colidx] = getReal(facets1->matrix[rowidx][colidx+1]);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons1_%d", rowidx);

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, nvars, vars, coefs,
            -getReal(facets1->matrix[rowidx][0]), SCIPinfinity(subscip)) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }
   for (rowidx = 0; rowidx < facets2->rowsize; ++rowidx)
   {
      for (colidx = 0; colidx < nvars; ++colidx)
         coefs[colidx] = getReal(facets2->matrix[rowidx][colidx+1]);

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "cons2_%d", rowidx);

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &cons, name, nvars, vars, coefs,
            -getReal(facets2->matrix[rowidx][0]), SCIPinfinity(subscip)) );
      SCIP_CALL( SCIPaddCons(subscip, cons) );
      SCIP_CALL( SCIPreleaseCons(subscip, &cons) );
   }

   SCIP_CALL( SCIPsolve(subscip) );

   if ( SCIPgetStatus(subscip) == SCIP_STATUS_OPTIMAL )
      intersect = TRUE;

   /* release variables */
   for (colidx = 0; colidx < nvars; ++colidx)
   {
      SCIP_CALL( SCIPreleaseVar(subscip, &vars[colidx]) );
   }

   SCIPfreeBlockMemoryArray(subscip, &coefs, nvars);
   SCIPfreeBlockMemoryArray(subscip, &vars, nvars);

   SCIP_CALL( SCIPfreeTransform(subscip) );
   SCIP_CALL( SCIPfree(&subscip) );

   return intersect;
}


/** checks whether a set of selected points from Y intersect with conv(X) */
static
SCIP_Bool intersects(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of points X */
   Datapoints*           Y,                  /**< set of points Y */
   int                   nX,                 /**< number of points in X */
   int                   nY,                 /**< number of points in Y */
   int*                  points,             /**< selected points */
   int                   npoints             /**< number of selected points */
   )
{
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhullY;
   dd_MatrixPtr facetsconvexhullX;
   SCIP_Bool success;
   SCIP_Bool intersects = FALSE;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( points != NULL );
   assert( npoints > 0 );

   /* compute convex hull of points */

   /* set up matrix of generators (points) */
   dd_set_global_constants();
   generators = constructGeneratorMatrixPoints(Y, points, npoints);
   facetsconvexhullY = computeConvexHullFacets(scip, generators, &success);

   dd_FreeMatrix(generators);

   if ( ! success )
   {
      dd_free_global_constants();
      return FALSE;
   }

   /* ignore convex sets that are not full-dimensional (lineality space is non-empty) */
   if ( set_card(facetsconvexhullY->linset) > 0 )
   {
      dd_FreeMatrix(facetsconvexhullY);
      dd_free_global_constants();
      return FALSE;
   }

   /* set up matrix of generators (X) */
   dd_set_global_constants();
   generators = constructGeneratorMatrixPoints(X, NULL, 0);
   facetsconvexhullX = computeConvexHullFacets(scip, generators, &success);

   dd_FreeMatrix(generators);

   if ( ! success )
   {
      dd_FreeMatrix(facetsconvexhullY);
      dd_free_global_constants();
      return FALSE;
   }

   /* check whether there exists a point in the intersection of both convex hull */
   intersects = polyhedraIntersect(facetsconvexhullX, facetsconvexhullY);

   dd_FreeMatrix(facetsconvexhullX);
   dd_FreeMatrix(facetsconvexhullY);
   dd_free_global_constants();

   return intersects;
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
   Datapoints* X;
   Datapoints* Y;
   int nX;
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

   assert( scip != NULL );
   assert( propdata != NULL );
   assert( nfixings != NULL );
   assert( infeasible != NULL );

   *nfixings = 0;
   *infeasible = FALSE;

   X = propdata->X;
   Y = propdata->Y;
   nX = propdata->nX;
   nY = propdata->nY;
   vars = propdata->vars;

   assert( X != NULL );
   assert( Y != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( vars != NULL );

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

   /* iterate over remaining infeasible points and check if adding them to separatedpoints
    * would intersect the convex hull of feasible points -> fix the variable to 0 in this case
    */
   for (i = 0; i < nY; ++i)
   {
      var = vars[i][idx];

      /* skip variables fixed to 0 */
      if ( SCIPvarGetUbLocal(var) < 0.5 )
         continue;

      /* fix already treated variables */
      if ( SCIPvarGetLbLocal(var) > 0.5 )
         continue;

      separatedpoints[nseparatedpoints] = i;

      if ( intersects(scip, X, Y, nX, nY, separatedpoints, nseparatedpoints + 1) )
      {
         SCIP_CALL( SCIPchgVarUb(scip, var, 0.0) );
         (*nfixings)++;
      }
   }

   SCIPfreeBufferArray(scip, &separatedpoints);

   return SCIP_OKAY;
}

/*
 * Callback methods of propagator
 */


/** initialization method of propagator (called after problem was transformed) */
static
SCIP_DECL_PROPINIT(propInitIntersection)
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

   SCIP_CALL( SCIPgetBoolParam(scip, "propagating/intersection/enabled", &propdata->enabled) );

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
SCIP_DECL_PROPFREE(propFreeIntersection)
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
SCIP_DECL_PROPPRESOL(propPresolIntersection)
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
      *result = SCIP_SUCCESS;
      *nfixedvars += nfixings;
   }

   return SCIP_OKAY;
}


/** execution method of propagator */
static
SCIP_DECL_PROPEXEC(propExecIntersection)
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

   if ( SCIPinProbing(scip) )
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
SCIP_DECL_PROPRESPROP(propRespropIntersection)
{  /*lint --e{715}*/
   assert( result != NULL );

   *result = SCIP_DIDNOTFIND;

   return SCIP_OKAY;
}


/*
 * propagator specific interface methods
 */

/** creates the intersection propagator and includes it in SCIP */
SCIP_RETCODE SCIPincludePropIntersection(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_PROP* prop;
   SCIP_PROPDATA* propdata;

   SCIP_CALL( propdataCreate(scip, &propdata) );

   /* include propagator */
   SCIP_CALL( SCIPincludePropBasic(scip, &prop, PROP_NAME, PROP_DESC, PROP_PRIORITY, PROP_FREQ, PROP_DELAY, PROP_TIMING,
         propExecIntersection, propdata) );

   assert( prop != NULL );

   /* set optional callbacks via setter functions */
   SCIP_CALL( SCIPsetPropInit(scip, prop, propInitIntersection) );
   SCIP_CALL( SCIPsetPropFree(scip, prop, propFreeIntersection) );
   SCIP_CALL( SCIPsetPropPresol(scip, prop, propPresolIntersection, PROP_PRESOL_PRIORITY,
         PROP_PRESOL_MAXROUNDS, PROP_PRESOLTIMING) );
   SCIP_CALL( SCIPsetPropResprop(scip, prop, propRespropIntersection) );

   SCIP_CALL( SCIPaddBoolParam(scip, "propagating/intersection/enabled",
         "whether the proagator is enabled", NULL, FALSE, TRUE, NULL, NULL) );

   return SCIP_OKAY;
}
