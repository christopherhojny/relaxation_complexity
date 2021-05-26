/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/*                                                                           */
/*    This file is part of the program colorbitopt                           */
/*                                                                           */
/*    an implementation of a branch-and-cut algorithm to solve the           */
/*    coloring problem by symmetry breaking methods based on orbitopes.      */
/*                                                                           */
/*    Copyright (C) 2005-2014  Marc Pfetsch                                  */
/*                                                                           */
/*                                                                           */
/*    Based on SCIP  --- Solving Constraint Integer Programs                 */
/*                                                                           */
/*    Copyright (C) 2002-2014 Konrad-Zuse-Zentrum                            */
/*                            fuer Informationstechnik Berlin                */
/*       mailto: scip@zib.de                                                 */
/*       SCIP is distributed under the terms of the SCIP Academic Licence,   */
/*       see file COPYING in the SCIP distribution.                          */
/*                                                                           */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file   getProblemName.cpp
 * @brief  Implementation of function to get the problem name from a filename
 * @author Marc Pfetsch
 */

#include "getProblemName.h"
#include <cassert>


/** Strip directory name and all suffixes and return result as string
 *
 *  This function is not portable and is likely to only work for UNIX systems. It, however, provides
 *  functionality that is not given by the function @c basename.
 */
std::string getProblemName(
   const char*           filename            //!< file name
   )
{
   // find last part of filename (used as problem name)
   unsigned int i = 0;
   while ( filename[i] != 0) // find end of string
      ++i;
   unsigned int l = i;       // end of string

   while ((i > 0) && (filename[i] != '.') && (filename[i] != '/')) // find last part
      --i;
   assert( i > 0 );

   if (filename[i] == '.')
   {
      l = i;
      while ((i > 0) && (filename[i] != '/'))   // find last part
	 --i;
   }

   if (filename[i] == '/')
      ++i;

   std::string problem_name;
   problem_name.reserve((unsigned int)(l + 1U));
   while ( (i < l) && (filename[i] != 0) )
      problem_name.push_back(filename[i++]);

   return problem_name;
}
