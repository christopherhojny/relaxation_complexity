/**@file   cons_conflict.h
 * @brief  Constraint handler to handle conflicts
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#ifndef __SCIP_CONS_CONFLICT_H__
#define __SCIP_CONS_CONFLICT_H__

#include "datapoints.h"

#ifdef __cplusplus
extern "C" {
#endif

/** creates the handler for conflict constraints and includes it in SCIP */
SCIP_RETCODE SCIPincludeConshdlrConflict(
   SCIP*                 scip                /**< SCIP data structure */
   );

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
   );

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
   );

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
   );

/** frees separation LP */
SCIP_RETCODE SCIPfreeSeparationLP(
   SCIP*                 subscip,            /**< subSCIP instance to be freed */
   int                   nfeaspoints,        /**< number of feasible points */
   int                   ninfpoints,         /**< number of selected infeasible points */
   int                   dim,                /**< dimension of problem */
   SCIP_VAR**            ineqvars,           /**< array of variables */
   SCIP_CONS**           validconss,         /**< array of valid conss */
   SCIP_CONS**           sepaconss           /**< array of separation conss */
   );


#ifdef __cplusplus
}
#endif

#endif
