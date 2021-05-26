// -*- C++ -*-
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

/**@file   getProblemName.h
 * @brief  Function to get the problem name from a filename
 * @author Marc Pfetsch
 */

#ifndef GETPROBLEMNAME_H
#define GETPROBLEMNAME_H

#include <string>

/** Strip directory name and all suffixes and return result as string
 *
 *  This function is not portable and is likely to only work for UNIX systems. It, however, provides
 *  functionality that is not given by the function @c basename.
 */
extern std::string getProblemName(
   const char*           filename            //!< file name
   );

#endif
