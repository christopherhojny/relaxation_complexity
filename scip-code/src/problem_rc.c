#include <stdio.h>

#include <scip/scip.h>
#include "scip/cons_indicator.h"
#include "scip/cons_linear.h"

#include "rcParams.h"
#include "rcPlugins.h"

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "pricer_pattern.h"
#include "problem_rc.h"
#include "probdata_rc_cg.h"
#include "probdata_rc_compact.h"
#include "typedefs.h"
#include "vardata_binpacking.h"


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
      display = TRUE;
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
      display = TRUE;
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


/** creates initial model for computing RC */
SCIP_RETCODE SCIPcreateModel(
   SCIP*                 scip,               /**< SCIP data structure */
   Datapoints*           X,                  /**< pointer to data points of X */
   Datapoints*           Y,                  /**< pointer to data points of Y */
   int*                  ub,                 /**< pointer to upper bound on RC(X,Y) */
   int                   lb,                 /**< lower bound on RC(X,Y) */
   int                   absmaxX,            /**< maximum absolute value of a coordinate in X */
   SCIP_Bool             secondphase         /**< whether we are in the second phase of the hybrid approach */
   )
{
   int method;

   assert( scip != NULL );
   assert( X != NULL );
   assert( Y != NULL );
   assert( ub != NULL );
   assert( *ub > 0 );

   SCIP_CALL( SCIPgetIntParam(scip, "rc/method", &method) );

   if ( method == METHOD_COMPACT_MIP || secondphase )
   {
      SCIP_CALL( SCIPprobdataCreateCompact(scip, "name", X, Y, ub, lb, absmaxX) );
   }
   else
   {
      SCIP_CALL( SCIPincludePricerPattern(scip) );
      SCIP_CALL( SCIPincludeBranchruleRyanFoster(scip) );
      SCIP_CALL( SCIPincludeConshdlrSamediff(scip) );
      SCIP_CALL( SCIPsetBoolParam(scip, "pricing/delvars", TRUE) );

      /* we have to ensure that SCIP does not remove the singleton sets, otherwise, the problem
       * might become infeasible after branching
       */
      if ( method == METHOD_CG )
      {
         SCIP_CALL( SCIPsetIntParam(scip, "constraints/setppc/maxprerounds", 0) );
      }

      SCIP_CALL( SCIPprobdataCreateCG(scip, "name", X, Y, absmaxX) );
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

   if ( method == METHOD_HYBRID )
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

      SCIP_VAR*** lhsvars;
      SCIP_VAR** rhsvars;
      SCIP_VAR*** violatedvars;
      SCIP_VAR** isusedvars;
      SCIP_CONS*** indicatorconss;
      SCIP_Real* inequality;
      SCIP_VAR* var;
      SCIP_Real slackval;
      SCIP_Real newtimelimit;
      int* pattern;
      int dimension;
      int len;

      /* solve the root node relaxation */
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", 1) );
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
      SCIP_CALL( SCIPsetRealParam(newscip, "limits/time", newtimelimit) );
      SCIP_CALL( SCIPsetLongintParam(newscip, "limits/nodes", -1) );

      newlb = (int) MAX(SCIPround(scip, dualbound), lb);
      newub = (int) MIN(SCIPround(scip, primalbound), ub);

      SCIP_CALL( SCIPcreateModel(newscip, X, Y, &newub, newlb, absmaxX, TRUE) );

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
               slackval = 0.001;

               solvar = SCIPprobdataGetVars(oldprobdata)[usedinequalities[j]];
               vardata = SCIPvarGetData(solvar);
               inequality = SCIPvardataGetInequality(vardata);

               /* compute -ay + b: note, stored inequality is already in shape b - ax >= 0 */
               for (k = 0; k < dimension; ++k)
                  slackval += inequality[k+1] * Y->points[i][k];
               slackval += inequality[0];

               SCIP_CALL( SCIPsetSolVal(newscip, newsol, var, MAX(slackval, 0.0)) );
            }
         }

         SCIPfreeBufferArray(scip, &usedinequalities);

         SCIP_CALL( SCIPaddSol(newscip, newsol, &stored) );
         SCIP_CALL( SCIPfreeSol(newscip, &newsol) );
      }

      SCIP_CALL( SCIPsolve(newscip) );
      SCIP_CALL( SCIPprintStatistics(newscip, NULL) );

      SCIP_CALL( writeSolutionCompact(newscip, filename, SCIPround(newscip, SCIPgetPrimalbound(newscip)), X->dimension) );

      /* free new structures */
      SCIP_CALL( SCIPfreeTransform(newscip) );
      SCIP_CALL( SCIPfree(&newscip) );
   }
   else
   {
      SCIP_CALL( SCIPsolve(scip) );

      /* output statistics */
      SCIPinfoMessage(scip, NULL, "\n");
      SCIP_CALL( SCIPprintStatistics(scip, NULL) );

      if ( method == METHOD_CG )
      {
         SCIP_CALL( writeSolutionCG(scip, filename) );
      }
      else
      {
         SCIP_CALL( writeSolutionCompact(scip, filename, ub, X->dimension) );
      }

   }

   return SCIP_OKAY;
}
