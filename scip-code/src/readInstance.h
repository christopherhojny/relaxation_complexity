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

/**@file   readInstance.h
 * @brief  function to read problem instance
 * @author Christopher Hojny
 */

// -*- C++ -*-

#ifndef READ_INSTANCE_H
#define READ_INSTANCE_H

#include <string>
#include "datapoints.h"

using namespace std;

/** reads a list of integer points from a file */
extern
SCIP_RETCODE readInstance(
   SCIP*                 scip,               /**< SCIP data structure */
   string                filename,           /**< name of file encoding data points */
   Datapoints**          datapoints,         /**< pointer to store data points */
   int*                  absmax,             /**< pointer to store maximum absolute value of coordinate (or NULL) */
   SCIP_Bool             printinfo           /**< Should information on instance be printed on screen? */
   );

#endif
