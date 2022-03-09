/**@file   branch_ryanfoster.c
 * @brief  Ryan/Foster branching rule
 * @author Christopher Hojny
 *
 * This file implements a Ryan/Foster-like branching rule.
 *
 * For branching, we are looking for two patterns P1 and P2 whose corresponding variables have fractional
 * values in the LP relaxation such that the intersection and symmetric difference of P1 and P2 are
 * non-empty. Then, we create two children by enforcing that P1 and P2 get the same or a different
 * value in the child nodes, respectively.
 */


#include <assert.h>
#include <string.h>

#include "branch_ryanfoster.h"
#include "cons_samediff.h"
#include "datapoints.h"
#include "probdata_rc_cg.h"
#include "vardata_binpacking.h"

#define BRANCHRULE_NAME            "RyanFoster"
#define BRANCHRULE_DESC            "Ryan/Foster branching rule"
#define BRANCHRULE_PRIORITY        50000
#define BRANCHRULE_MAXDEPTH        -1
#define BRANCHRULE_MAXBOUNDDIST    1.0

/** branching execution method for fractional LP solutions */
static
SCIP_DECL_BRANCHEXECLP(branchExeclpRyanFoster)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   Datapoints* Y;
   int nY;
   SCIP_VAR** lpcands;
   SCIP_Real* lpcandsfrac;
   int nlpcands;

   SCIP_NODE* childsame;
   SCIP_NODE* childdiffer;
   SCIP_CONS* conssame;
   SCIP_CONS* consdiffer;

   SCIP_CONSHDLR* samediffconshdlr;
   SCIP_CONS** samediffconss;
   int nsamediffconss;
   int* branchingpairs;
   int npairs;
   int dimension;

   int id1 = -1;
   int id2 = -1;

   int i;
   int j;
   int v;

   SCIP_Real maxval = -1.0;
   int maxid1;
   int maxid2;

   assert( scip != NULL );
   assert( branchrule != NULL );
   assert( strcmp(SCIPbranchruleGetName(branchrule), BRANCHRULE_NAME) == 0 );
   assert( result != NULL );

   *result = SCIP_DIDNOTRUN;

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   Y = SCIPprobdataGetY(probdata);
   assert( Y != NULL );

   nY = Y->ndatapoints;
   dimension = Y->dimension;

   /* get fractional LP candidates */
   SCIP_CALL( SCIPgetLPBranchCands(scip, &lpcands, NULL, &lpcandsfrac, NULL, &nlpcands, NULL) );
   assert(nlpcands > 0);

   /* replace fractionality by fractionality w.r.t. 0.5 */
   for (i = 0; i < nlpcands; ++i)
      lpcandsfrac[i] = 0.5 - MIN(lpcandsfrac[i], 1 - lpcandsfrac[i]);

   /* sort candidates w.r.t. fractionality */
   SCIPsortRealPtr(lpcandsfrac, (void**) lpcands, nlpcands);

   samediffconshdlr = SCIPfindConshdlr(scip, "samediff");
   assert( samediffconshdlr != NULL );

   samediffconss = SCIPconshdlrGetConss(samediffconshdlr);
   nsamediffconss = SCIPconshdlrGetNConss(samediffconshdlr);

   SCIP_CALL( SCIPallocBufferArray(scip, &branchingpairs, 2 * nsamediffconss) );
   npairs = 0;

   for (i = 0; i < nsamediffconss; ++i)
   {
      SCIP_CONS* cons;

      cons = samediffconss[i];

      /* ignore inactive constraints (not laying on the path from the current node to the root) */
      if( !SCIPconsIsActive(cons) )
         continue;

      branchingpairs[npairs++] = SCIPgetItemid1Samediff(scip, cons);
      branchingpairs[npairs++] = SCIPgetItemid2Samediff(scip, cons);
   }
   npairs = npairs / 2;

   /* find pair of sets with maximum common value in LP solution */
   for (i = 0; i < nlpcands; ++i)
   {
      SCIP_Real* inequ1;

      inequ1 = SCIPvardataGetInequality(SCIPvarGetData(lpcands[i]));

      for (j = i + 1; j < nlpcands; ++j)
      {
         SCIP_Real* inequ2;

         inequ2 = SCIPvardataGetInequality(SCIPvarGetData(lpcands[j]));

         /* we can stop early if we cannot exceed the best known value maxval */
         if ( SCIPisLE(scip, lpcandsfrac[i] + lpcandsfrac[j], maxval) )
            break;

         /* check whether sets i and j have non-empty intersection and symmetric difference */
         for (v = 0; v < nY; ++v)
         {
            SCIP_Real val1a;
            SCIP_Real val1b;
            int k;
            int w;

            /* evaluate inequality */
            val1a = inequ1[0];
            for (k = 0; k < dimension; ++k)
               val1a += inequ1[k+1] * Y->points[v][k];

            for (w = v + 1; w < nY; ++w)
            {
               SCIP_Real val2a;
               SCIP_Real val2b;

               /* skip already handled pairs */
               for (k = 0; k < npairs; ++k)
               {
                  if ( (branchingpairs[2*k] == v && branchingpairs[2*k+1] == w) ||
                     (branchingpairs[2*k+1] == v && branchingpairs[2*k] == w) )
                     break;
               }

               if ( k < npairs )
                  continue;

               /* evaluate inequality */
               val1b = inequ1[0];
               for (k = 0; k < dimension; ++k)
                  val1b += inequ1[k+1] * Y->points[w][k];
               val2a = inequ2[0];
               for (k = 0; k < dimension; ++k)
                  val2a += inequ2[k+1] * Y->points[v][k];
               val2b = inequ2[0];
               for (k = 0; k < dimension; ++k)
                  val2b += inequ2[k+1] * Y->points[w][k];

               /* find elements in intersection (id1) and symmetric difference (id2) (if they exist) */
               if ( SCIPisLT(scip, val1a, 0.0) && SCIPisLT(scip, val2a, 0.0) )
                  id1 = v;
               if ( SCIPisLT(scip, val1b, 0.0) && SCIPisLT(scip, val2b, 0.0) )
                  id1 = w;
               if ( SCIPisGE(scip, val1a, 0.0) && SCIPisLT(scip, val2a, 0.0) )
                  id2 = v;
               if ( SCIPisLT(scip, val1a, 0.0) && SCIPisGE(scip, val2a, 0.0) )
                  id2 = v;
               if ( SCIPisGE(scip, val1b, 0.0) && SCIPisLT(scip, val2b, 0.0) )
                  id2 = w;
               if ( SCIPisLT(scip, val1b, 0.0) && SCIPisGE(scip, val2b, 0.0) )
                  id2 = w;

               if ( id1 != -1 && id2 != -1 )
               {
                  /* keep track of elements in intersection and symmetric difference is sets
                   * yield larger fractional value
                   */
                  if ( SCIPisLT(scip, maxval, lpcandsfrac[i] + lpcandsfrac[j]) )
                  {
                     maxval = lpcandsfrac[i] + lpcandsfrac[j];
                     maxid1 = id1;
                     maxid2 = id2;
                     id1 = -1;
                     id2 = -1;
                  }

                  break;
               }

               id1 = -1;
               id2 = -1;
            }

            if ( w < nY )
               break;
         }

         if ( v < nY )
            break;
      }
   }

   /* create the branch-and-bound tree child nodes of the current node */
   SCIP_CALL( SCIPcreateChild(scip, &childsame, 0.0, SCIPgetLocalTransEstimate(scip)) );
   SCIP_CALL( SCIPcreateChild(scip, &childdiffer, 0.0, SCIPgetLocalTransEstimate(scip)) );

   /* create corresponding constraints */
   SCIP_CALL( SCIPcreateConsSamediff(scip, &conssame, "same", maxid1, maxid2, SAME, childsame, TRUE) );
   SCIP_CALL( SCIPcreateConsSamediff(scip, &consdiffer, "differ", maxid1, maxid2, DIFFER, childdiffer, TRUE) );

  /* add constraints to nodes */
   SCIP_CALL( SCIPaddConsNode(scip, childsame, conssame, NULL) );
   SCIP_CALL( SCIPaddConsNode(scip, childdiffer, consdiffer, NULL) );

   /* release constraints */
   SCIP_CALL( SCIPreleaseCons(scip, &conssame) );
   SCIP_CALL( SCIPreleaseCons(scip, &consdiffer) );

   SCIPfreeBufferArray(scip, &branchingpairs);

  *result = SCIP_BRANCHED;

   return SCIP_OKAY;
}

/** creates the Ryan/Foster-like branching rule and includes it in SCIP */
SCIP_RETCODE SCIPincludeBranchruleRyanFoster(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_BRANCHRULEDATA* branchruledata;
   SCIP_BRANCHRULE* branchrule;

   /* create ryan foster branching rule data */
   branchruledata = NULL;
   branchrule = NULL;

   /* include branching rule */
   SCIP_CALL( SCIPincludeBranchruleBasic(scip, &branchrule, BRANCHRULE_NAME, BRANCHRULE_DESC,
         BRANCHRULE_PRIORITY, BRANCHRULE_MAXDEPTH, BRANCHRULE_MAXBOUNDDIST, branchruledata) );

   assert( branchrule != NULL );

   SCIP_CALL( SCIPsetBranchruleExecLp(scip, branchrule, branchExeclpRyanFoster) );

   return SCIP_OKAY;
}
