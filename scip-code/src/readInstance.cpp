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

/**@file   readInstance.cpp
 * @brief  function to read problem instance
 * @author Christopher Hojny
 */

#include <iostream>
#include <fstream>
#include <sstream>

#include "datapoints.h"
#include "readInstance.h"

using namespace std;

/** reads a list of integer points from a file */
SCIP_RETCODE readInstance(
   SCIP*                 scip,               /**< SCIP data structure */
   string                filename,           /**< name of file encoding data points */
   Datapoints**          datapoints,         /**< pointer to store data points */
   int*                  absmax,             /**< pointer to store maximum absolute value of coordinate (or NULL) */
   SCIP_Bool             printinfo           /**< Should information on instance be printed on screen? */
   )
{
   int ndatapoints;
   int dimension;
   int val;

   assert( scip != NULL );
   assert( datapoints != NULL );

   cout << "Reading data points from file" << endl;
   ifstream input(filename);
   if (!input)
   {
      cout << "No input file found." << endl;
   }

   SCIP_CALL( SCIPallocBlockMemory(scip, datapoints) );

   // read number of data points and dimension of data points
   input >> ndatapoints;
   input >> dimension;

   (*datapoints)->ndatapoints = ndatapoints;
   (*datapoints)->dimension = dimension;

   if ( absmax != NULL )
      *absmax = 0;

   cout << "Read file containing " << ndatapoints << " data points of dimension " << dimension << endl;

   // allocate memory for data points
   SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->points, ndatapoints) );
   for (int i = 0; i < ndatapoints; ++i)
   {
      SCIP_CALL( SCIPallocBlockMemoryArray(scip, &(*datapoints)->points[i], dimension) );
   }

   // read information from file
   for (int i = 0; i < ndatapoints; ++i)
   {
      for (int j = 0; j < dimension; ++j)
      {
         input >> val;
         (*datapoints)->points[i][j] = val;

         if ( absmax != NULL )
         {
            if ( ABS(val) > *absmax )
               *absmax = ABS(val);
         }

         if ( printinfo )
            printf(" %d", val);
      }
      if ( printinfo )
         printf("\n");
   }

   return SCIP_OKAY;;
}
