/**@file   vardata_compact.c
 * @author Christopher Hojny
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/

#include "vardata_compact.h"

/** Variable data which is attached to all variables indicating whether an inequality is violated */
struct SCIP_VarData
{
   int                   idx;                /**< index of corresponding inequality */
};

/**@name Local methods
 *
 * @{
 */

/** create a vardata */
static
SCIP_RETCODE vardataCreate(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int                   idx                 /**< index of corresponding inequality */
   )
{
   SCIP_CALL( SCIPallocBlockMemory(scip, vardata) );

   (*vardata)->idx = idx;

   return SCIP_OKAY;
}

/** frees user data of variable */
static
SCIP_RETCODE vardataDelete(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata             /**< vardata to delete */
   )
{
   SCIPfreeBlockMemory(scip, vardata);

   return SCIP_OKAY;
}

/**@} */


/**@name Callback methods
 *
 * @{
 */

/** frees user data of transformed variable (called when the transformed variable is freed) */
static
SCIP_DECL_VARDELTRANS(vardataDelTrans)
{
   SCIP_CALL( vardataDelete(scip, vardata) );

   return SCIP_OKAY;
}/*lint !e715*/

/**@} */


/**@name Interface methods
 *
 * @{
 */

/** creates is violated variable */
SCIP_RETCODE SCIPcreateVarIsviolatedCompact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VAR**            var,                /**< pointer to variable object */
   const char*           name,               /**< name of variable, or NULL for automatic name creation */
   SCIP_Real             obj,                /**< objective function value */
   SCIP_VARDATA*         vardata             /**< user data for this specific variable */
   )
{
   assert(scip != NULL);
   assert(var != NULL);

   /* create a basic variable object */
   SCIP_CALL( SCIPcreateVarBasic(scip, var, name, 0.0, 1.0, obj, SCIP_VARTYPE_BINARY) );
   assert(*var != NULL);

   /* set callback functions */
   SCIPvarSetData(*var, vardata);
   SCIPvarSetDeltransData(*var, vardataDelTrans);

   return SCIP_OKAY;
}

/** create variable data */
SCIP_RETCODE SCIPvardataCreateCompact(
   SCIP*                 scip,               /**< SCIP data structure */
   SCIP_VARDATA**        vardata,            /**< pointer to vardata */
   int                   idx                 /**< index of corresponding inequality */
   )
{
   SCIP_CALL( vardataCreate(scip, vardata, idx) );

   return SCIP_OKAY;
}

/** get index of corresponding inequality */
int SCIPvardataGetInequalityidx(
   SCIP_VARDATA*         vardata             /**< variable data */
   )
{
   return vardata->idx;
}

/**@} */
