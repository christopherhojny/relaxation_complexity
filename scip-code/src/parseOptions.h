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

/**@file   parse_options.h
 * @brief  parser for options
 * @author Christopher Hojny
 */

#ifndef PARSEOPTIONS_H
#define PARSEOPTIONS_H

#include <vector>
#include <string>
#include <cstring>


/** reads command line options */
bool readOptions(
   int                   argc,               //!< number of arguments
   const char*const*     argv,               //!< array of arguments
   int                   nminargs,           //!< minimum number of necessary arguments
   int                   nmaxargs,           //!< maximum number of necessary arguments
   std::vector<std::string>& mainargs,       //!< necessary arguments
   std::vector<std::string>& options,        //!< options enabled by (-<option letter>)
   std::vector<std::string> knownOptions     //!< known options (option letters)
   )
{
   assert( nminargs <= nmaxargs );

   if ( argc <= nminargs )
   {
      printf("Expected at least %d and at most %d arguments, but got %d.\n", nminargs, nmaxargs, argc - 1);
      return false;
   }

   mainargs.clear();
   options.clear();

   // read necessary arguments
   for (int i = 1; i <= nminargs; ++i)
   {
      if ( argv[i][0] == '-' )
      {
         printf("Expected necessary argument, but got %s.\n", argv[i]);
         return false;
      }
      mainargs.push_back(std::string(argv[i]));
   }

   int i = nminargs + 1;
   for (i = nminargs + 1; i <= MIN(nmaxargs, argc-1); ++i)
   {
      // found optional argument, break
      if ( argv[i][0] == '-' )
         break;
      mainargs.push_back(std::string(argv[i]));
   }

   // read options
   for ( ; i < argc; ++i)
   {
      if ( argv[i][0] != '-' )
      {
         printf("expected optional argument, but got %s.\n", argv[i]);
         return false;
      }

      std::string arg;
      arg.push_back(argv[i][1]);

      bool found = false;
      for (long unsigned int j = 0; j < knownOptions.size(); ++j)
      {
         if ( knownOptions[j] == arg )
         {
            found = true;
            break;
         }
      }
      if ( ! found )
      {
         printf("Got unknown option -%s.\n", arg.c_str());
         return false;
      }

      options.push_back(arg);
      ++i;

      if ( i == argc )
      {
         printf("expected value of optional argument %s, but got nothing.\n", argv[i-1]);
         return false;
      }

      options.push_back(argv[i]);
   }

   return true;
}

#endif
