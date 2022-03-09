#include "compute_symmetries.h"

#include <unordered_map>

/* include bliss graph */
#include <bliss/defs.hh>
#include <bliss/graph.hh>

/** struct for bliss callback */
struct BLISS_Data
{
   SCIP*                 scip;               /**< SCIP pointer */
   int                   npoints;            /**< number of points permutations are acting on */
   int                   nperms;             /**< number of permutations */
   int**                 perms;              /**< permutation generators as (nperms x npermvars) matrix */
   int                   nmaxperms;          /**< maximal number of permutations */
};


/** callback function for bliss */
static
void blisshook(
   void*                 user_param,         /**< parameter supplied at call to bliss */
   unsigned int          n,                  /**< size of aut vector */
   const unsigned int*   aut                 /**< automorphism */
   )
{
   assert( aut != NULL );
   assert( user_param != NULL );

   BLISS_Data* data = (BLISS_Data*) user_param;
   assert( data->scip != NULL );
   assert( data->npoints < (int) n );

   /* copy first part of automorphism */
   bool isIdentity = true;
   int* p = 0;
   if ( SCIPallocBlockMemoryArray(data->scip, &p, data->npoints) != SCIP_OKAY )
      return;

   for (int j = 0; j < data->npoints; ++j)
   {
      /* convert index of variable-level 0-nodes to variable indices */
      p[j] = (int) aut[j];
      if ( p[j] != j )
         isIdentity = false;
   }

   /* ignore trivial generators, i.e. generators that only permute the constraints */
   if ( isIdentity )
   {
      SCIPfreeBlockMemoryArray(data->scip, &p, data->npoints);
      return;
   }

   /* check whether we should allocate space for perms */
   if ( data->nmaxperms <= 0 )
   {
      data->nmaxperms = 100;   /* seems to cover many cases */

      if ( SCIPallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms) != SCIP_OKAY )
         return;
   }
   else if ( data->nperms >= data->nmaxperms )    /* check whether we need to resize */
   {
      int newsize = SCIPcalcMemGrowSize(data->scip, data->nperms + 1);
      assert( newsize >= data->nperms );

      if ( SCIPreallocBlockMemoryArray(data->scip, &data->perms, data->nmaxperms, newsize) != SCIP_OKAY )
         return;

      data->nmaxperms = newsize;
   }

   data->perms[data->nperms++] = p;
}


/** gets information on data points for a specific coordinate */
static
SCIP_RETCODE getInfo(
   SCIP_Real**           Xpoints,            /**< points in point class X */
   SCIP_Real**           Ypoints,            /**< points in point class Y */
   int                   nX,                 /**< number of points in X */
   int                   nY,                 /**< number of points in Y */
   int                   coord,              /**< coordinate of points to be checked */
   int*                  minval,             /**< minimum value of coord in X and Y */
   int*                  maxval              /**< maximum value of coord in X and Y */
   )
{
   assert( Xpoints != NULL );
   assert( Ypoints != NULL );
   assert( nX > 0 );
   assert( nY > 0 );
   assert( minval != NULL );
   assert( maxval != NULL );

   *minval = (int) Ypoints[0][coord];
   *maxval = (int) Ypoints[0][coord];
   for (int i = 0; i < nY; ++i)
   {
      if ( Ypoints[i][coord] > *maxval )
         *maxval = Ypoints[i][coord];
      else if ( Ypoints[i][coord] < *minval )
         *minval = Ypoints[i][coord];
   }
   for (int i = 0; i < nX; ++i)
   {
      if ( Xpoints[i][coord] > *maxval )
         *maxval = Xpoints[i][coord];
      else if ( Xpoints[i][coord] < *minval )
         *minval = Xpoints[i][coord];
   }

   return SCIP_OKAY;
}


/** computes a key value for a coordinate in a certain dimension */
int getKeyCoordinate(
   int                   coordinateval,      /**< coordinatevalue of a points */
   int                   minval,             /**< minimum coordinate value */
   int                   maxval,             /**< maximum coordinate value */
   int                   dim                 /**< dimension of point */
   )
{
   return coordinateval - minval + dim * maxval;
}


/** computes common symmetries of X and Y */
SCIP_RETCODE computeSymmetries(
   SCIP*                 scip,               /**< SCIP instance */
   Datapoints*           X,                  /**< set of feasible points */
   Datapoints*           Y,                  /**< set of infeasible points */
   int***                perms,              /**< pointer to store permutations */
   int*                  nperms,             /**< pointer to store number of permutations stored in perms */
   int*                  nmaxperms           /**< pointer to store maximum number of permutations fiiting in perms */
   )
{
   SCIP_Real** Xpoints;
   SCIP_Real** Ypoints;
   int nX;
   int nY;
   int dim;
   int minval;
   int maxval;

   Xpoints = X->points;
   Ypoints = Y->points;
   nX = X->ndatapoints;
   nY = Y->ndatapoints;
   dim = X->dimension;

   /* We create a colored complete bipartite graph with partition P, C, where P encodes the points and C the coordinates.
    * Points from Y are colored with INT_MAX, points from X with INT_MAX - 1, and, for each possible coordinate value
    * p_i of a point in P, we create a coordinate node with the corresponding color. Instead of absolute coordinates,
    * we use relative coordinates to compute symmetries up to translations of the point set. For example, if coordinate i
    * attains values between -3 and 1, we create nodes for coordinate i that are colored 0, ..., 4.
    * There is an edge between p in P and c_i(colored k) if p_i = k.
    *
    * Automorphisms of this graph restricted to Y correspond to reorderings of the points in Y that leave the
    * problem invariant.
    **/

   /* create bliss graph */
   bliss::Graph G(0);

   /* create a node for every point in Y with color INT_MAX */
   for (int v = 0; v < nY; ++v)
      (void) G.add_vertex(INT_MAX);

   /* create a node for every point in X with color INT_MAX - 1 */
   for (int v = 0; v < nY; ++v)
      (void) G.add_vertex(INT_MAX - 1);

   for (int j = 0; j < dim; ++j)
   {
      /* compute coordinate nodes and connect points with these nodes */
      SCIP_CALL( getInfo(Xpoints, Ypoints, nX, nY, j, &minval, &maxval) );

      int* coordinates;
      SCIP_CALL( SCIPallocBufferArray(scip, &coordinates, maxval - minval + 1) );

      /* init used array */
      for (int i = 0; i < maxval - minval + 1; ++i)
         coordinates[i] = -1;

      /* for each point in X/Y, create node for coordinate j with color according to coordinate value */
      for (int i = 0; i < nX; ++i)
      {
         int val = Xpoints[i][j] - minval;

         if ( coordinates[val] < 0 )
            coordinates[val] = (int) G.add_vertex(val);

         G.add_edge(nY + i, coordinates[val]);
      }
      for (int i = 0; i < nY; ++i)
      {
         int val = Ypoints[i][j] - minval;

         if ( coordinates[val] < 0 )
            coordinates[val] = (int) G.add_vertex(val);

         G.add_edge(i, coordinates[val]);
      }

      SCIPfreeBufferArray(scip, &coordinates);
   }

   /* compute automorphisms */
   bliss::Stats stats;
   BLISS_Data data;
   data.scip = scip;
   data.npoints = nY;
   data.nperms = 0;
   data.nmaxperms = 0;
   data.perms = NULL;

   /* Prefer splitting partition cells corresponding to variables over those corresponding
    * to inequalities. This is because we are only interested in the action
    * of the automorphism group on the variables, and we need a base for this action */
   G.set_splitting_heuristic(bliss::Graph::shs_f);
   /* disable component recursion as advised by Tommi Junttila from bliss */
   G.set_component_recursion(false);

   /* start search */
   G.find_automorphisms(stats, blisshook, (void*) &data);

   if ( data.nperms > 0 )
   {
      *nperms = data.nperms;
      *perms = data.perms;
      *nmaxperms = data.nmaxperms;
   }
   else
   {
      *nperms = 0;
      *nmaxperms = 0;
      *perms = NULL;
   }

   return SCIP_OKAY;
}

