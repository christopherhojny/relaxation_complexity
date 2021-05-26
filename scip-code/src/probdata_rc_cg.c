/**@file   probdata_rc_compact.c
 * @brief  Problem data for computing RC using a compact model
 * @author Christopher Hojny
 *
 * This file handles the main problem data used in that project.
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <string.h>

#include "datapoints.h"
#include "pricer_pattern.h"
#include "probdata_rc_cg.h"
#include "vardata_binpacking.h"

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "auxiliary_cdd.h"
#include "math.h"
#include "convex_hull.h"
#include "typedefs.h"

#include "scip/scip.h"
#include "scip/scipdefplugins.h"
#include "scip/cons_linear.h"
#include "scip/cons_indicator.h"
#include "scip/cons_setppc.h"

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
   int                   dimension;          /**< dimension of data points */
   int                   absmaxX;            /**< maximum absolute entry of a coordinate in X */

   /* data of the model */
   SCIP_VAR**            modelvars;          /**< variables related to sets in the covering formulation */
   int                   nmodelvars;         /**< number of model variables */
   int                   maxnmodelvars;      /**< maximum number of variables that can be stored in modelvars */
   SCIP_CONS**           coverconss;         /**< covering constraints */
};


/**@name Event handler properties
 *
 * @{
 */

#define EVENTHDLR_NAME         "addedvar"
#define EVENTHDLR_DESC         "event handler for catching added variables"

/**@} */

/**@name Callback methods of event handler
 *
 * @{
 */

/** execution method of event handler */
static
SCIP_DECL_EVENTEXEC(eventExecAddedVar)
{  /*lint --e{715}*/
   assert(eventhdlr != NULL);
   assert(strcmp(SCIPeventhdlrGetName(eventhdlr), EVENTHDLR_NAME) == 0);
   assert(event != NULL);
   assert(SCIPeventGetType(event) == SCIP_EVENTTYPE_VARADDED);

   SCIPdebugMsg(scip, "exec method of event handler for added variable to probdata\n");

   /* add new variable to probdata */
   SCIP_CALL( SCIPprobdataAddVar(scip, SCIPgetProbData(scip), SCIPeventGetVar(event)) );

   return SCIP_OKAY;
}

/**@} */

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
   int                   absmaxX,            /**< maximum absolute entry of a coordinate in X */
   SCIP_VAR**            vars,               /**< variables related to sets in the covering formulation */
   int                   nvars   ,           /**< number of variables in vars */
   int                   maxnvars,           /**< maximum number of variables that can be stored in modelvars */
   SCIP_CONS**           coverconss          /**< covering constraints */
   )
{
   int nX;
   int nY;
   int dimension;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( X != NULL );
   assert( Y != NULL );

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
   (*probdata)->absmaxX = absmaxX;
   (*probdata)->nmodelvars = nvars;
   (*probdata)->maxnmodelvars = maxnvars;

   /* possible copy variable arrays */
   if ( vars != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->modelvars, vars, maxnvars) );
   }
   else
      (*probdata)->modelvars = NULL;

   /* duplicate constraint arrays */
   if ( coverconss != NULL )
   {
      SCIP_CALL( SCIPduplicateBlockMemoryArray(scip, &(*probdata)->coverconss, coverconss, nY) );
   }
   else
      (*probdata)->coverconss = NULL;

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
   int nY;

   assert( scip != NULL );
   assert( probdata != NULL );
   assert( (*probdata)->X != NULL );
   assert( (*probdata)->Y != NULL );

   nY = (*probdata)->nY;

   /* release all variables */
   for (i = 0; i < (*probdata)->nmodelvars; ++i)
   {
      SCIP_CALL( SCIPreleaseVar(scip, &(*probdata)->modelvars[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->modelvars, (*probdata)->maxnmodelvars);

   /* release all constraints */
   for (i = 0; i < nY; ++i)
   {
      SCIP_CALL( SCIPreleaseCons(scip, &(*probdata)->coverconss[i]) );
   }
   SCIPfreeBlockMemoryArrayNull(scip, &(*probdata)->coverconss, nY);

   /* free probdata */
   SCIPfreeBlockMemory(scip, probdata);

   return SCIP_OKAY;
}


/** finds inequality separating a given point from a set of points */
static
SCIP_RETCODE findSeparatingInequality(
   Datapoints*           X,                  /**< set of points */
   int                   absmaxX,            /**< maximum absolute value of coordinate in X */
   SCIP_Real*            point,              /**< point to be separated */
   int                   dimension,          /**< dimension of point */
   SCIP_Real*            inequality,         /**< allocated memory to store inequality */
   SCIP_Real             epsilon             /**< the epsilon used to rate an inequality as violated */
   )
{
   char name[SCIP_MAXSTRLEN];
   SCIP* subscip;
   SCIP_VAR** inequvars;
   SCIP_Real* tmpvals;
   SCIP_CONS** validconss;
   SCIP_CONS* sepacons;
   int nX;
   int i;
   int j;

   assert( X != NULL );
   assert( point != NULL );
   assert( dimension > 0 );
   assert( X->dimension == dimension );
   assert( inequality != NULL );

   nX = X->ndatapoints;

   /* basic setup */
   SCIP_CALL( SCIPcreate(&subscip) );
   SCIP_CALL( SCIPcreateProbBasic(subscip, "findInequality") );

   /* do not abort subproblem on CTRL-C */
   SCIP_CALL( SCIPsetBoolParam(subscip, "misc/catchctrlc", FALSE) );

   /* disable output to console */
   SCIP_CALL( SCIPsetIntParam(subscip, "display/verblevel", 0) );

   /* load plugins */
   SCIP_CALL( SCIPincludeDefaultPlugins(subscip) );

   /* create variables */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &inequvars, dimension + 1) );
   for (i = 0; i < dimension; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "a_%d", i);

      SCIP_CALL( SCIPcreateVar(subscip, &inequvars[i], name, -1.0, 1.0, 0.0,
            SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
      SCIP_CALL( SCIPaddVar(subscip, inequvars[i]) );
   }
   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "b");

   SCIP_CALL( SCIPcreateVar(subscip, &inequvars[dimension], name,
         -dimension * absmaxX, dimension * absmaxX, 0.0,
         SCIP_VARTYPE_CONTINUOUS, TRUE, FALSE, NULL, NULL, NULL, NULL, NULL) );
   SCIP_CALL( SCIPaddVar(subscip, inequvars[dimension]) );

   /* create constraints */
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &tmpvals, dimension + 1 ) );
   SCIP_CALL( SCIPallocBlockMemoryArray(subscip, &validconss, nX) );
   tmpvals[dimension] = -1.0;
   for (i = 0; i < nX; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "validcons_%d", i);

      for (j = 0; j < dimension; ++j)
         tmpvals[j] = X->points[i][j];

      SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &validconss[i], name,
            dimension + 1, inequvars, tmpvals, -SCIPinfinity(subscip), 0.0) );
      SCIP_CALL( SCIPaddCons(subscip, validconss[i]) );
      SCIP_CALL( SCIPreleaseCons(subscip, &validconss[i]) );
   }

   (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "sepacons_%d", i);

   for (j = 0; j < dimension; ++j)
      tmpvals[j] = -point[j];
   tmpvals[dimension] = 1.0;

   SCIP_CALL( SCIPcreateConsBasicLinear(subscip, &sepacons, name,
         dimension + 1,  inequvars, tmpvals, -SCIPinfinity(subscip), -epsilon) );
   SCIP_CALL( SCIPaddCons(subscip, sepacons) );
   SCIP_CALL( SCIPreleaseCons(subscip, &sepacons) );

   SCIP_CALL( SCIPsolve(subscip) );

   inequality[0] = SCIPgetSolVal(subscip, SCIPgetBestSol(subscip), inequvars[dimension]);
   for (i = 0; i < dimension; ++i)
      inequality[i+1] = -SCIPgetSolVal(subscip, SCIPgetBestSol(subscip), inequvars[i]);

   for (i = 0; i <= dimension; ++i)
      SCIP_CALL( SCIPreleaseVar(subscip, &inequvars[i]) );

   SCIPfreeBlockMemoryArray(subscip, &validconss, nX);
   SCIPfreeBlockMemoryArray(subscip, &tmpvals, dimension + 1);
   SCIPfreeBlockMemoryArray(subscip, &inequvars, dimension + 1);

   SCIP_CALL( SCIPfreeTransform(subscip) );
   SCIP_CALL( SCIPfree(&subscip) );


   return SCIP_OKAY;
}


/** creates the initial variables of the problem */
static
SCIP_RETCODE SCIPcreateInitialVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata            /**< problem data */
   )
{
   char name[SCIP_MAXSTRLEN];
   dd_MatrixPtr generators;
   dd_MatrixPtr facetsconvexhull;
   dd_rowrange rowidx;
   dd_colrange colidx;
   SCIP_VAR* newvar;
   SCIP_Real* inequality;
   SCIP_Bool success;
   SCIP_Real eps;
   int nY;
   int* set;
   int lenset;
   int dimension;
   int method;
   int i;

   assert( scip != NULL );
   assert( probdata != NULL );

   assert( probdata->nX > 0 );
   assert( probdata->nY > 0 );
   assert( probdata->dimension > 0 );

   nY = probdata->nY;
   dimension = probdata->dimension;

   SCIP_CALL( SCIPgetRealParam(scip, "rc/epsilon", &eps) );

   /* compute facets of convex hull of X*/
   dd_set_global_constants();
   generators = constructGeneratorMatrixPoints(probdata->X, NULL, 0);

   facetsconvexhull = computeConvexHullFacets(scip, generators, &success);

   dd_FreeMatrix(generators);

   if ( ! success )
   {
      dd_free_global_constants();

      return SCIP_ERROR;
   }

   /* compute the sets of points in Y separated by the facet defining inequalities */
   SCIP_CALL( SCIPallocBufferArray(scip, &set, probdata->Y->ndatapoints) );
   SCIP_CALL( SCIPallocBufferArray(scip, &inequality, dimension + 1) );

   for (rowidx = 0; rowidx < facetsconvexhull->rowsize; ++rowidx)
   {
      SCIP_VARDATA* vardata;
      SCIP_Real absmax = 1.0;

      lenset = 0;

      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "init_%d", rowidx);

      for (i = 0; i < nY; ++i)
      {
         SCIP_Real val;

         /* evaluate the i-th point in facedt rowidx */
         val = getReal(facetsconvexhull->matrix[rowidx][0]);
         for (colidx = 1; colidx <= dimension; ++colidx)
            val += getReal(facetsconvexhull->matrix[rowidx][colidx]) * probdata->Y->points[i][colidx - 1];

         /* the point is separated by the inequality, add it to the current set */
         if ( SCIPisLT(scip, val, 0.0) )
            set[lenset++] = i;
      }

      /* create variable for the packing pattern corresponding to facet */
      SCIP_CALL( SCIPcreateVarBinpacking(scip, &newvar, name, 1.0, TRUE, TRUE, TRUE, NULL) );

      /* add variable to the problem */
      SCIP_CALL( SCIPaddVar(scip, newvar) );

      /* store variable in the problme data */
      SCIP_CALL( SCIPprobdataAddVar(scip, probdata, newvar) );

      /* add variable to corresponding set covering constraints */
      for (i = 0; i < lenset; ++i)
      {
         SCIP_CALL( SCIPaddCoefSetppc(scip, probdata->coverconss[set[i]], newvar) );
      }

      /* create the variable data for the variable; the variable data contains the information in which constraints the
       * variable appears */
      for (i = 0; i <= dimension; ++i)
      {
	  if ( ABS(getReal(facetsconvexhull->matrix[rowidx][i])) > absmax )
	     absmax = ABS(getReal(facetsconvexhull->matrix[rowidx][i]));
      }
      for (i = 0; i <= dimension; ++i)
         inequality[i] = getReal(facetsconvexhull->matrix[rowidx][i]) / absmax;
      SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, set, lenset, inequality, dimension + 1) );

      /* add the variable data to the variable */
      SCIPvarSetData(newvar, vardata);

      SCIP_CALL( SCIPchgVarUbLazy(scip, newvar, 1.0) );

      /* release variable */
      SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
   }

   /* create singleton sets to avoid infeasible problems after branching */
   SCIP_CALL( SCIPgetIntParam(scip, "rc/method", &method) );
   if ( method == METHOD_CG )
   {
      for (i = 0; i < nY; ++i)
      {
         SCIP_VARDATA* vardata;

         (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "singleton_%d", i);

         /* create variable for the packing pattern which contains only this item */
         SCIP_CALL( SCIPcreateVarBinpacking(scip, &newvar, name, 1.0, TRUE, FALSE, FALSE, NULL) );

         /* add variable to the problem */
         SCIP_CALL( SCIPaddVar(scip, newvar) );

         /* store variable in the problme data */
         SCIP_CALL( SCIPprobdataAddVar(scip, probdata, newvar) );

         /* add variable to corresponding set covering constraints */
         SCIP_CALL( SCIPaddCoefSetppc(scip, probdata->coverconss[i], newvar) );

         set[0] = i;

         /* add a dummy inequality */
         SCIP_CALL( findSeparatingInequality(probdata->X, probdata->absmaxX, probdata->Y->points[i],
               dimension, inequality, eps) );

         SCIP_CALL( SCIPvardataCreateBinpacking(scip, &vardata, set, 1, inequality, dimension + 1) );

         /* add the variable data to the variable */
         SCIPvarSetData(newvar, vardata);

         SCIP_CALL( SCIPchgVarUbLazy(scip, newvar, 1.0) );

         /* release variable */
         SCIP_CALL( SCIPreleaseVar(scip, &newvar) );
      }
   }

   SCIPfreeBufferArray(scip, &inequality);
   SCIPfreeBufferArray(scip, &set);

   dd_FreeMatrix(facetsconvexhull);
   dd_free_global_constants();

   return SCIP_OKAY;
}


/**@} */

/**@name Callback methods of problem data
 *
 * @{
 */

/** frees user data of original problem (called when the original problem is freed) */
static
SCIP_DECL_PROBDELORIG(probdelorigRCCG)
{
   SCIPdebugMsg(scip, "free original problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** creates user data of transformed problem by transforming the original user problem data
 *  (called after problem was transformed) */
static
SCIP_DECL_PROBTRANS(probtransRCCG)
{
   int nY;
   int nvars;

   nY = sourcedata->nY;
   nvars = sourcedata->nmodelvars;

   /* create transform probdata */
   SCIP_CALL( probdataCreate(scip, targetdata, sourcedata->X, sourcedata->Y, sourcedata->absmaxX,
         sourcedata->modelvars, sourcedata->nmodelvars, sourcedata->maxnmodelvars,
         sourcedata->coverconss) );

   /* transform all constraints */
   SCIP_CALL( SCIPtransformConss(scip, nY, (*targetdata)->coverconss, (*targetdata)->coverconss) );

   /* transform all variables */
   SCIP_CALL( SCIPtransformVars(scip, nvars, (*targetdata)->modelvars, (*targetdata)->modelvars) );

   return SCIP_OKAY;
}

/** frees user data of transformed problem (called when the transformed problem is freed) */
static
SCIP_DECL_PROBDELTRANS(probdeltransRCCG)
{
   SCIPdebugMsg(scip, "free transformed problem data\n");

   SCIP_CALL( probdataFree(scip, probdata) );

   return SCIP_OKAY;
}

/** solving process initialization method of transformed data (called before the branch and bound process begins) */
static
SCIP_DECL_PROBINITSOL(probinitsolRCCG)
{
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* catch variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPcatchEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, NULL) );

   return SCIP_OKAY;
}

/** solving process deinitialization method of transformed data (called before the branch and bound data is freed) */
static
SCIP_DECL_PROBEXITSOL(probexitsolRCCG)
{  /*lint --e{715}*/
   SCIP_EVENTHDLR* eventhdlr;

   assert(probdata != NULL);

   /* drop variable added event */
   eventhdlr = SCIPfindEventhdlr(scip, "addedvar");
   assert(eventhdlr != NULL);

   SCIP_CALL( SCIPdropEvent(scip, SCIP_EVENTTYPE_VARADDED, eventhdlr, NULL, -1) );


   return SCIP_OKAY;
}

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** sets up the problem data */
SCIP_RETCODE SCIPprobdataCreateCG(
   SCIP*                 scip,               /**< SCIP data structure */
   const char*           probname,           /**< problem name */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int                   absmaxX             /**< maximum absolute value of a coordinate in X */
   )
{
   SCIP_PROBDATA* probdata;
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONS** conss;
   char name[SCIP_MAXSTRLEN];
   int i;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( absmaxX > 0 );

   /* create event handler if it does not exist yet */
   if( SCIPfindEventhdlr(scip, EVENTHDLR_NAME) == NULL )
   {
      SCIP_CALL( SCIPincludeEventhdlrBasic(scip, NULL, EVENTHDLR_NAME, EVENTHDLR_DESC, eventExecAddedVar, NULL) );
   }

   /* create problem in SCIP and add non-NULL callbacks via setter functions */
   SCIP_CALL( SCIPcreateProbBasic(scip, probname) );

   SCIP_CALL( SCIPsetProbDelorig(scip, probdelorigRCCG) );
   SCIP_CALL( SCIPsetProbTrans(scip, probtransRCCG) );
   SCIP_CALL( SCIPsetProbDeltrans(scip, probdeltransRCCG) );
   SCIP_CALL( SCIPsetProbInitsol(scip, probinitsolRCCG) );
   SCIP_CALL( SCIPsetProbExitsol(scip, probexitsolRCCG) );

   /* set objective sense */
   SCIP_CALL( SCIPsetObjsense(scip, SCIP_OBJSENSE_MINIMIZE) );
   SCIP_CALL( SCIPsetObjIntegral(scip) );

   SCIP_CALL( SCIPallocBufferArray(scip, &conss, Y->ndatapoints) );

   /* create set covering constraints for each item */
   for (i = 0; i < Y->ndatapoints; ++i)
   {
      (void) SCIPsnprintf(name, SCIP_MAXSTRLEN, "covercons_%d", i);

      SCIP_CALL( SCIPcreateConsBasicSetcover(scip, &conss[i], name, 0, NULL) );

      /* declare constraint modifiable for adding variables during pricing */
      SCIP_CALL( SCIPsetConsModifiable(scip, conss[i], TRUE) );
      SCIP_CALL( SCIPaddCons(scip, conss[i]) );
   }

   /* create problem data */
   SCIP_CALL( probdataCreate(scip, &probdata, X, Y, absmaxX,
         NULL, 0, 0, conss) );

   SCIP_CALL( SCIPcreateInitialVariables(scip, probdata) );

   /* set user problem data */
   SCIP_CALL( SCIPsetProbData(scip, probdata) );

   conshdlr = SCIPfindConshdlr(scip, "samediff");
   SCIP_CALL( SCIPpricerPatternActivate(scip, X, Y, absmaxX, conss, Y->ndatapoints, conshdlr) );

   /* free local buffer arrays */
   SCIPfreeBufferArray(scip, &conss);

   return SCIP_OKAY;
}


/** adds given variable to the problem data */
SCIP_RETCODE SCIPprobdataAddVar(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_PROBDATA*        probdata,           /**< problem data */
   SCIP_VAR*             var                 /**< variables to add */
   )
{
   /* check if enough memory is left */
   if( probdata->maxnmodelvars == probdata->nmodelvars )
   {
      int newsize;
      newsize = MAX(100, probdata->maxnmodelvars * 2);
      SCIP_CALL( SCIPreallocBlockMemoryArray(scip, &probdata->modelvars, probdata->maxnmodelvars, newsize) );
      probdata->maxnmodelvars = newsize;
   }

   /* caputure variables */
   SCIP_CALL( SCIPcaptureVar(scip, var) );

   probdata->modelvars[probdata->nmodelvars] = var;
   probdata->nmodelvars++;

   return SCIP_OKAY;
}


/** returns covering constraints of problem */
SCIP_CONS** SCIPprobdataCGGetCoverconss(
   SCIP_PROBDATA*        probdata            /** problem data */
   )
{
   assert( probdata != NULL );

   return probdata->coverconss;
}


/** returns number of covering constraints in problem */
int SCIPprobdataCGGetNCoverconss(
   SCIP_PROBDATA*        probdata            /** problem data */
  )
{
   assert( probdata != NULL );

   return probdata->nY;
}


/** returns number of variables in problem */
int SCIPprobdataGetNVars(
   SCIP_PROBDATA*        probdata            /** problem data */
   )
{
   assert( probdata != NULL );

   return probdata->nmodelvars;
}


/** returns variables in problem */
SCIP_VAR** SCIPprobdataGetVars(
   SCIP_PROBDATA*        probdata            /** problem data */
   )
{
   assert( probdata != NULL );

   return probdata->modelvars;
}

/** returns variables in problem */
Datapoints* SCIPprobdataGetY(
   SCIP_PROBDATA*        probdata            /** problem data */
   )
{
   assert( probdata != NULL );

   return probdata->Y;
}

/**@} */
