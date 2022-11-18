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

/**@file   auxiliary_cdd.c
 * @brief  file implementing auxiliary functions for CDD
 * @author Christopher Hojny
 */

#include "cddlib/setoper.h"
#include "cddlib/cddmp.h"
#include "cddlib/cdd.h"
#include "math.h"


#include "auxiliary_cdd.h"

/*  modified copy of dd_WriteReal() from cddio.c */
double getReal(
   mytype x
   )
{
   long ix1;
   long ix2;
   long ix;
   double ax;

   ax=dd_get_d(x);
   ix1= (long) (fabs(ax) * 10000. + 0.5);
   ix2= (long) (fabs(ax) + 0.5);
   ix2= ix2*10000;
   if ( ix1 == ix2 )
   {
      if ( dd_Positive(x) )
      {
         ix = (long)(ax + 0.5);
      }
      else
      {
         ix = (long)(-ax + 0.5);
         ix = -ix;
      }
      return ix;
   }

  return ax;
}
