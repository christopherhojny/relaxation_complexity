/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*                  This file is part of the program and library             */
/*         SCIP --- Solving Constraint Integer Programs                      */
/*                                                                           */
/*  Copyright 2002-2022 Zuse Institute Berlin                                */
/*                                                                           */
/*  Licensed under the Apache License, Version 2.0 (the "License");          */
/*  you may not use this file except in compliance with the License.         */
/*  You may obtain a copy of the License at                                  */
/*                                                                           */
/*      http://www.apache.org/licenses/LICENSE-2.0                           */
/*                                                                           */
/*  Unless required by applicable law or agreed to in writing, software      */
/*  distributed under the License is distributed on an "AS IS" BASIS,        */
/*  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. */
/*  See the License for the specific language governing permissions and      */
/*  limitations under the License.                                           */
/*                                                                           */
/*  You should have received a copy of the Apache-2.0 license                */
/*  along with SCIP; see the file LICENSE. If not visit scipopt.org.         */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**
 * CHANGES TO ORIGINAL FILE:
 *
 * May 20, 2021
 * - removed definition and calls of consdataPrint() and consdataCheck()
 * - removed implementation of CONSPRINT method
 */

/**@file   cons_samediff.c
 * @brief  Constraint handler stores the local branching decision data
 * @author Timo Berthold
 * @author Stefan Heinz
 *
 * This constraint handler is used to store the branching decision of the \ref BINPACKING_BRANCHING "Ryan/Foster branching rule"
 * which is implemented in \ref branch_ryanfoster.c.
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include <assert.h>
#include <string.h>

#include "cons_samediff.h"
#include "probdata_rc_cg.h"
#include "vardata_binpacking.h"


/**@name Constraint handler properties
 *
 * @{
 */

#define CONSHDLR_NAME          "samediff"
#define CONSHDLR_DESC          "stores the local branching decisions"
#define CONSHDLR_ENFOPRIORITY         0 /**< priority of the constraint handler for constraint enforcing */
#define CONSHDLR_CHECKPRIORITY  9999999 /**< priority of the constraint handler for checking feasibility */
#define CONSHDLR_PROPFREQ             1 /**< frequency for propagating domains; zero means only preprocessing propagation */
#define CONSHDLR_EAGERFREQ            1 /**< frequency for using all instead of only the useful constraints in separation,
                                         *   propagation and enforcement, -1 for no eager evaluations, 0 for first only */
#define CONSHDLR_DELAYPROP        FALSE /**< should propagation method be delayed, if other propagators found reductions? */
#define CONSHDLR_NEEDSCONS         TRUE /**< should the constraint handler be skipped, if no constraints are available? */

#define CONSHDLR_PROP_TIMING       SCIP_PROPTIMING_BEFORELP

/**@} */

/*
 * Data structures
 */

/** Constraint data for  \ref cons_samediff.c "SameDiff" constraints */
struct SCIP_ConsData
{
   int                   itemid1;            /**< item id one */
   int                   itemid2;            /**< item id two */
   CONSTYPE              type;               /**< stores whether the items have to be in the SAME or DIFFER packing */
   int                   npropagatedvars;    /**< number of variables that existed, the last time, the related node was
                                              *   propagated, used to determine whether the constraint should be
                                              *   repropagated*/
   int                   npropagations;      /**< stores the number propagations runs of this constraint */
   unsigned int          propagated:1;       /**< is constraint already propagated? */
   SCIP_NODE*            node;               /**< the node in the B&B-tree at which the cons is sticking */
};

/**@name Local methods
 *
 * @{
 */

/** create constraint data */
static
SCIP_RETCODE consdataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata,           /**< pointer to store the constraint data */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node                /**< the node in the B&B-tree at which the cons is sticking */
   )
{
   assert( scip != NULL );
   assert( consdata != NULL );
   assert( itemid1 >= 0 );
   assert( itemid2 >= 0 );
   assert( type == DIFFER || type == SAME );

   SCIP_CALL( SCIPallocBlockMemory(scip, consdata) );

   (*consdata)->itemid1 = itemid1;
   (*consdata)->itemid2 = itemid2;
   (*consdata)->type = type;
   (*consdata)->npropagatedvars = 0;
   (*consdata)->npropagations = 0;
   (*consdata)->propagated = FALSE;
   (*consdata)->node = node;

   return SCIP_OKAY;
}


/** fixes a variable to zero if the corresponding packings are not valid for this constraint/node (due to branching) */
static
SCIP_RETCODE checkVariable(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR*             var,                /**< variables to check  */
   int*                  nfixedvars,         /**< pointer to store the number of fixed variables */
   SCIP_Bool*            cutoff              /**< pointer to store if a cutoff was detected */
   )
{
   SCIP_VARDATA* vardata;
   int* consids;
   int nconsids;

   SCIP_Bool existid1;
   SCIP_Bool existid2;
   CONSTYPE type;

   SCIP_Bool fixed;
   SCIP_Bool infeasible;

   int pos;

   assert(scip != NULL);
   assert(consdata != NULL);
   assert(var != NULL);
   assert(nfixedvars != NULL);
   assert(cutoff != NULL);

   /* if variables is locally fixed to zero continue */
   if( SCIPvarGetUbLocal(var) < 0.5 )
      return SCIP_OKAY;

   /* check if the packing which corresponds to the variable is feasible for this constraint */
   vardata = SCIPvarGetData(var);

   nconsids = SCIPvardataGetNConsids(vardata);
   consids = SCIPvardataGetConsids(vardata);

   existid1 = SCIPsortedvecFindInt(consids, consdata->itemid1, nconsids, &pos);
   existid2 = SCIPsortedvecFindInt(consids, consdata->itemid2, nconsids, &pos);
   type = consdata->type;

   if( (type == SAME && existid1 != existid2) || (type == DIFFER && existid1 && existid2) )
   {
      SCIP_CALL( SCIPfixVar(scip, var, 0.0, &infeasible, &fixed) );

      if( infeasible )
      {
         assert( SCIPvarGetLbLocal(var) > 0.5 );
         SCIPdebugMsg(scip, "-> cutoff\n");
         (*cutoff) = TRUE;
      }
      else
      {
         assert(fixed);
         (*nfixedvars)++;
      }
   }

   return SCIP_OKAY;
}

/** fixes variables to zero if the corresponding packings are not valid for this sonstraint/node (due to branching) */
static
SCIP_RETCODE consdataFixVariables(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA*        consdata,           /**< constraint data */
   SCIP_VAR**            vars,               /**< generated variables */
   int                   nvars,              /**< number of generated variables */
   SCIP_RESULT*          result              /**< pointer to store the result of the fixing */
   )
{
   int nfixedvars;
   int v;
   SCIP_Bool cutoff;

   nfixedvars = 0;
   cutoff = FALSE;

   SCIPdebugMsg(scip, "check variables %d to %d\n", consdata->npropagatedvars, nvars);

   for( v = consdata->npropagatedvars; v < nvars && !cutoff; ++v )
   {
      SCIP_CALL( checkVariable(scip, consdata, vars[v], &nfixedvars, &cutoff) );
   }

   SCIPdebugMsg(scip, "fixed %d variables locally\n", nfixedvars);

   if( cutoff )
      *result = SCIP_CUTOFF;
   else if( nfixedvars > 0 )
      *result = SCIP_REDUCEDDOM;

   return SCIP_OKAY;
}


/** frees samediff constraint data */
static
SCIP_RETCODE consdataFree(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONSDATA**       consdata            /**< pointer to the constraint data */
   )
{
   assert(consdata != NULL);
   assert(*consdata != NULL);

   SCIPfreeBlockMemory(scip, consdata);

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees specific constraint data */
static
SCIP_DECL_CONSDELETE(consDeleteSamediff)
{  /*lint --e{715}*/
   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(consdata != NULL);
   assert(*consdata != NULL);

   /* free samediff constraint */
   SCIP_CALL( consdataFree(scip, consdata) );

   return SCIP_OKAY;
}

/** transforms constraint data into data belonging to the transformed problem */
static
SCIP_DECL_CONSTRANS(consTransSamediff)
{  /*lint --e{715}*/
   SCIP_CONSDATA* sourcedata;
   SCIP_CONSDATA* targetdata;

   assert(conshdlr != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(SCIPgetStage(scip) == SCIP_STAGE_TRANSFORMING);
   assert(sourcecons != NULL);
   assert(targetcons != NULL);

   sourcedata = SCIPconsGetData(sourcecons);
   assert(sourcedata != NULL);

   /* create constraint data for target constraint */
   SCIP_CALL( consdataCreate(scip, &targetdata,
         sourcedata->itemid1, sourcedata->itemid2, sourcedata->type, sourcedata->node) );

   /* create target constraint */
   SCIP_CALL( SCIPcreateCons(scip, targetcons, SCIPconsGetName(sourcecons), conshdlr, targetdata,
         SCIPconsIsInitial(sourcecons), SCIPconsIsSeparated(sourcecons), SCIPconsIsEnforced(sourcecons),
         SCIPconsIsChecked(sourcecons), SCIPconsIsPropagated(sourcecons),
         SCIPconsIsLocal(sourcecons), SCIPconsIsModifiable(sourcecons),
         SCIPconsIsDynamic(sourcecons), SCIPconsIsRemovable(sourcecons), SCIPconsIsStickingAtNode(sourcecons)) );

   return SCIP_OKAY;
}

/** constraint enforcing method of constraint handler for LP solutions */
#define consEnfolpSamediff NULL

/** constraint enforcing method of constraint handler for pseudo solutions */
#define consEnfopsSamediff NULL

/** feasibility check method of constraint handler for integral solutions */
#define consCheckSamediff NULL

/** domain propagation method of constraint handler */
static
SCIP_DECL_CONSPROP(consPropSamediff)
{  /*lint --e{715}*/
   SCIP_PROBDATA* probdata;
   SCIP_CONSDATA* consdata;

   SCIP_VAR** vars;
   int nvars;
   int c;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(result != NULL);

   SCIPdebugMsg(scip, "propagation constraints of constraint handler <"CONSHDLR_NAME">\n");

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   vars = SCIPprobdataGetVars(probdata);
   nvars = SCIPprobdataGetNVars(probdata);

   *result = SCIP_DIDNOTFIND;

   for( c = 0; c < nconss; ++c )
   {
      consdata = SCIPconsGetData(conss[c]);

#ifndef NDEBUG
      {
         /* check if there are no equal consdatas */
         SCIP_CONSDATA* consdata2;
         int i;

         for( i = c+1; i < nconss; ++i )
         {
            consdata2 = SCIPconsGetData(conss[i]);
            assert( !(consdata->itemid1 == consdata2->itemid1
                  && consdata->itemid2 == consdata2->itemid2
                  && consdata->type == consdata2->type) );
            assert( !(consdata->itemid1 == consdata2->itemid2
                  && consdata->itemid2 == consdata2->itemid1
                  && consdata->type == consdata2->type) );
         }
      }
#endif

      if( !consdata->propagated )
      {
         SCIPdebugMsg(scip, "propagate constraint <%s> ", SCIPconsGetName(conss[c]));
         SCIPdebug( consdataPrint(scip, consdata, NULL) );

         SCIP_CALL( consdataFixVariables(scip, consdata, vars, nvars, result) );
         consdata->npropagations++;

         if( *result != SCIP_CUTOFF )
         {
            consdata->propagated = TRUE;
            consdata->npropagatedvars = nvars;
         }
         else
            break;
      }
   }

   return SCIP_OKAY;
}

/** variable rounding lock method of constraint handler */
#define consLockSamediff NULL

/** constraint activation notification method of constraint handler */
static
SCIP_DECL_CONSACTIVE(consActiveSamediff)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->npropagatedvars <= SCIPprobdataGetNVars(probdata));

   SCIPdebugMsg(scip, "activate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
      SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
   SCIPdebug( consdataPrint(scip, consdata, NULL) );

   if( consdata->npropagatedvars != SCIPprobdataGetNVars(probdata) )
   {
      SCIPdebugMsg(scip, "-> mark constraint to be repropagated\n");
      consdata->propagated = FALSE;
      SCIP_CALL( SCIPrepropagateNode(scip, consdata->node) );
   }

   return SCIP_OKAY;
}

/** constraint deactivation notification method of constraint handler */
static
SCIP_DECL_CONSDEACTIVE(consDeactiveSamediff)
{  /*lint --e{715}*/
   SCIP_CONSDATA* consdata;
   SCIP_PROBDATA* probdata;

   assert(scip != NULL);
   assert(strcmp(SCIPconshdlrGetName(conshdlr), CONSHDLR_NAME) == 0);
   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);
   assert(consdata->propagated || SCIPgetNChildren(scip) == 0);

   probdata = SCIPgetProbData(scip);
   assert(probdata != NULL);

   SCIPdebugMsg(scip, "deactivate constraint <%s> at node <%"SCIP_LONGINT_FORMAT"> in depth <%d>: ",
      SCIPconsGetName(cons), SCIPnodeGetNumber(consdata->node), SCIPnodeGetDepth(consdata->node));
   SCIPdebug( consdataPrint(scip, consdata, NULL) );

   /* set the number of propagated variables to current number of variables is SCIP */
   consdata->npropagatedvars = SCIPprobdataGetNVars(probdata);

   return SCIP_OKAY;
}


/**@} */

/**@name Interface methods
 *
 * @{
 */

/** creates the handler for samediff constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrSamediff(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   SCIP_CONSHDLRDATA* conshdlrdata = NULL;
   SCIP_CONSHDLR* conshdlr = NULL;

   /* include constraint handler */
   SCIP_CALL( SCIPincludeConshdlrBasic(scip, &conshdlr, CONSHDLR_NAME, CONSHDLR_DESC,
         CONSHDLR_ENFOPRIORITY, CONSHDLR_CHECKPRIORITY, CONSHDLR_EAGERFREQ, CONSHDLR_NEEDSCONS,
         consEnfolpSamediff, consEnfopsSamediff, consCheckSamediff, consLockSamediff,
         conshdlrdata) );
   assert(conshdlr != NULL);

   SCIP_CALL( SCIPsetConshdlrDelete(scip, conshdlr, consDeleteSamediff) );
   SCIP_CALL( SCIPsetConshdlrTrans(scip, conshdlr, consTransSamediff) );
   SCIP_CALL( SCIPsetConshdlrProp(scip, conshdlr, consPropSamediff, CONSHDLR_PROPFREQ, CONSHDLR_DELAYPROP,
         CONSHDLR_PROP_TIMING) );
   SCIP_CALL( SCIPsetConshdlrActive(scip, conshdlr, consActiveSamediff) );
   SCIP_CALL( SCIPsetConshdlrDeactive(scip, conshdlr, consDeactiveSamediff) );

   return SCIP_OKAY;
}

/** creates and captures a samediff constraint */
SCIP_RETCODE SCIPcreateConsSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS**           cons,               /**< pointer to hold the created constraint */
   const char*           name,               /**< name of constraint */
   int                   itemid1,            /**< item id one */
   int                   itemid2,            /**< item id two */
   CONSTYPE              type,               /**< stores whether the items have to be in the SAME or DIFFER packing */
   SCIP_NODE*            node,               /**< the node in the B&B-tree at which the cons is sticking */
   SCIP_Bool             local               /**< is constraint only valid locally? */
   )
{
   SCIP_CONSHDLR* conshdlr;
   SCIP_CONSDATA* consdata;

   /* find the samediff constraint handler */
   conshdlr = SCIPfindConshdlr(scip, CONSHDLR_NAME);
   if( conshdlr == NULL )
   {
      SCIPerrorMessage("samediff constraint handler not found\n");
      return SCIP_PLUGINNOTFOUND;
   }

   /* create the constraint specific data */
   SCIP_CALL( consdataCreate(scip, &consdata, itemid1, itemid2, type, node) );

   /* create constraint */
   SCIP_CALL( SCIPcreateCons(scip, cons, name, conshdlr, consdata, FALSE, FALSE, FALSE, FALSE, TRUE,
         local, FALSE, FALSE, FALSE, TRUE) );

   SCIPdebugMsg(scip, "created constraint: ");
   SCIPdebug( consdataPrint(scip, consdata, NULL) );

   return SCIP_OKAY;
}

/** returns item id one */
int SCIPgetItemid1Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->itemid1;
}

/** returns item id two */
int SCIPgetItemid2Samediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->itemid2;
}

/** return constraint type SAME or DIFFER */
CONSTYPE SCIPgetTypeSamediff(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_CONS*            cons                /**< samediff constraint */
   )
{
   SCIP_CONSDATA* consdata;

   assert(cons != NULL);

   consdata = SCIPconsGetData(cons);
   assert(consdata != NULL);

   return consdata->type;
}

/**@} */
