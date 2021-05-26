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
