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

/**@file   rcParams.cpp
 * @brief  set up parameters for computing RC
 * @author Christopher Hojny
 */

#include "rcParams.h"


/* default settings for parameters - for documentation see below */

/* information regarding problem */
#define DEFAULT_PRINTDATAINFO                FALSE
#define DEFAULT_HANDLESYMMETRY                TRUE
#define DEFAULT_USEORBITOPE                   TRUE
#define DEFAULT_USESYMRESACKS                 TRUE
#define DEFAULT_HIDINGSETCUTS                 TRUE
#define DEFAULT_METHOD                           2
#define DEFAULT_EPSILON                      0.001
#define DEFAULT_CUTOFFEPSILON               0.0001
#define DEFAULT_MAXSOLPRICER                     1
#define DEFAULT_PREPROCESSING                FALSE
#define DEFAULT_MAXSEPA                      FALSE
#define DEFAULT_HEURNODELIMIT                 1000

/** Set basic SCIP parameters that are relevant for computing RC */
extern
SCIP_RETCODE setSCIPParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   return SCIP_OKAY;
}


/** Introduce parameters that are relevant for computing RC */
extern
SCIP_RETCODE addRCParameters(
   SCIP*                 scip                /**< SCIP data structure */
   )
{
   assert( scip != NULL );

   /* information regarding problem */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/printdatainfo", "Shall information on the data points be printed?",
         NULL, TRUE, DEFAULT_PRINTDATAINFO, NULL, NULL) );

   /* type of method to solve problem */
   SCIP_CALL( SCIPaddIntParam(scip, "rc/method", "Which method shall be used to find RC (0: compact MIP, 1: column generation, 2: hybrid CG/compact, 3: conflict MIP, 4: hybrid CG/conflict )",
         NULL, TRUE, DEFAULT_METHOD, 0, 4, NULL, NULL) );

   /* whether we shall use symmetry handling techniques */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/handlesymmetry", "Shall symmetry be handled in compact model?",
         NULL, TRUE, DEFAULT_HANDLESYMMETRY, NULL, NULL) );

   /* whether we shall use orbitopes when symmetry handling is active  */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/useorbitope", "Shall orbitopes be used if symmetry handling is active?",
         NULL, TRUE, DEFAULT_USEORBITOPE, NULL, NULL) );

   /* whether we shall use symresacks to handle Y symmetries when symmetry handling is active  */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/usesymresacks", "Shall symresacks be used if symmetry handling is active to handle Y-symmetries?",
         NULL, TRUE, DEFAULT_USESYMRESACKS, NULL, NULL) );

   /* whether we shall use cuts based on hiding sets */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/hidingsetcuts", "Shall hiding set based cuts be added to the model?",
         NULL, TRUE, DEFAULT_HIDINGSETCUTS, NULL, NULL) );

   /* the epsilon used to rate an inequality as violated */
   SCIP_CALL( SCIPaddRealParam(scip, "rc/epsilon", "The epsilon used to rate an inequality as violated",
         NULL, TRUE, DEFAULT_EPSILON, 0.00001, 1.0, NULL, NULL) );

   /* the epsilon used to rate an inequality as violated */
   SCIP_CALL( SCIPaddRealParam(scip, "rc/cutoffepsilon", "The epsilon used to accept solutions in the pricing problem.",
         NULL, TRUE, DEFAULT_CUTOFFEPSILON, 0.0, 1.0, NULL, NULL) );

   /* number of solutions excepted in the pricing subproblem */
   SCIP_CALL( SCIPaddIntParam(scip, "rc/maxsolpricer", "Maximum number of solutions above the cutoff limit after which pricing terminates early.",
         NULL, TRUE, DEFAULT_MAXSOLPRICER, 1, INT_MAX, NULL, NULL) );

   /* shall we perform preprocessing (remove nonObservers from Y) */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/preprocessing", "shall we perform preprocessing (remove nonObservers from Y)?",
         NULL, TRUE, DEFAULT_PREPROCESSING, NULL, NULL) );

   /* whether the pricer shall separate maximal sets only  */
   SCIP_CALL( SCIPaddBoolParam(scip, "rc/maxsepa", "whether the pricer shall separate maximal sets only",
         NULL, TRUE, DEFAULT_PREPROCESSING, NULL, NULL) );

   /* maximum number of B&B nodes used in heuristics per iteration */
   SCIP_CALL( SCIPaddLongintParam(scip, "rc/heurnodelimit", "maximum number of B&B nodes used in heuristics per iteration",
         NULL, TRUE, DEFAULT_HEURNODELIMIT, -1, SCIP_LONGINT_MAX, NULL, NULL) );

   /* maximum number of B&B nodes used in pricing problems to improve separated inequality */
   SCIP_CALL( SCIPaddLongintParam(scip, "rc/pricerimprovenodelimit",
         "maximum number of B&B nodes used in pricing problems to improve separated inequality",
         NULL, TRUE, DEFAULT_HEURNODELIMIT, -1, SCIP_LONGINT_MAX, NULL, NULL) );

   return SCIP_OKAY;
}
