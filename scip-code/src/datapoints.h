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

/**@file   datapoints.h
 * @brief Declaration of data points
 * @author Christopher Hojny
 *
 * This file contains the definitions of a struct to store data points
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef DATAPOINTS_H
#define DATAPOINTS_H

#include "scip/scip.h"

#ifdef __cplusplus
extern "C" {
#endif

/** data of a node */
typedef struct Datapoints
{
   int                   dimension;          /**< dimension of set of data points */
   int                   ndatapoints;        /**< number of data points */
   SCIP_Real**           points;             /**< (ndatapoints x dimension)-array of data points */
} Datapoints;

#ifdef __cplusplus
}
#endif

#endif

