/**@file   main.cpp
 * @brief  main file for computing the relaxation complexity
 * @author Christopher Hojny
 */

#include "rcParams.h"
#include "rcPlugins.h"
#include "datapoints.h"
#include "getProblemName.h"
#include "parseOptions.h"
#include "readInstance.h"
#include "problem_rc.h"
#include "probdata_rc_cg.h"
#include "probdata_rc_compact.h"
#include "vardata_binpacking.h"
#include "typedefs.h"

#include <scip/scip.h>

#include <iostream>
#include <fstream>
#include <limits>

/**@brief computeRC(X,Y) */
static
SCIP_RETCODE computeRC(
   std::string           filenameX,          //!< name of file encoding points X
   std::string           filenameY,          //!< name of file encoding points Y
   int                   ub,                 //!< upper bound on RC(X,Y)
   int                   lb,                 //!< lower bound on RC(X,Y)
   char*                 relaxfilename,      //!< filename to store relaxation
   const char*           settings,           //!< Possible name of setting file
   double                timeLimit,          //!< time limit
   double                memLimit,           //!< memory limit
   SCIP_Longint          nodeLimit,          //!< node limit
   int                   displayFreq         //!< display frequency
   )
{  /*lint --e{429}*/

   // initialize SCIP
   SCIP* scip;
   SCIP_CALL( SCIPcreate(&scip) );

   // output SCIP banner
#ifdef SCIP_THREADSAFE_MESSAGEHDLRS
   SCIPprintVersion(scip, NULL);
#else
   SCIPprintVersion(NULL);
#endif
   std::cout << "\n" << std::endl;   // (to force flush)

   // load basic plugins
   SCIP_CALL( includeRCPlugins(scip) );

   // add our own parameters
   SCIP_CALL( setSCIPParameters(scip) );
   SCIP_CALL( addRCParameters(scip) );

   if ( settings != 0 )
   {
      if ( ! SCIPfileExists(settings) )
      {
         SCIPerrorMessage("Setting file <%s> does not exist.\n\n", settings);
      }
      else
      {
         SCIPinfoMessage(scip, 0, "Reading parameters from <%s>.\n\n", settings);
         SCIP_CALL( SCIPreadParams(scip, settings) );
      }
   }

   // initialize problem class
   std::string problemName = getProblemName(filenameX.c_str());

   /* output changed parameters */
   SCIPinfoMessage(scip, 0, "Changed settings:\n");
   SCIP_CALL( SCIPwriteParams(scip, 0, FALSE, TRUE) );
   SCIPinfoMessage(scip, 0, "\n");

   // set limits
   if ( timeLimit < 1e20 )
   {
      SCIPinfoMessage(scip, 0, "Setting time limit to %g.\n", timeLimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/time", timeLimit) );
   }
   if ( memLimit < 1e20 )
   {
      SCIPinfoMessage(scip, 0, "Setting memory limit to %g.\n", memLimit);
      SCIP_CALL( SCIPsetRealParam(scip, "limits/memory", memLimit) );
   }
   if ( nodeLimit < SCIP_LONGINT_MAX )
   {
      SCIPinfoMessage(scip, 0, "Setting node limit to %lld.\n", nodeLimit);
      SCIP_CALL( SCIPsetLongintParam(scip, "limits/nodes", nodeLimit) );
   }
   if ( displayFreq < INT_MAX )
   {
      SCIPinfoMessage(scip, 0, "Setting display frequency to %d.\n", displayFreq);
      SCIP_CALL( SCIPsetIntParam(scip, "display/freq", displayFreq) );
   }
   SCIPinfoMessage(scip, 0, "\n");

   // create problem
   SCIP_Bool printinfo;
   Datapoints* X;
   Datapoints* Y;
   int absmaxX;
   SCIP_CALL( SCIPgetBoolParam(scip, "rc/printdatainfo", &printinfo) );
   SCIP_CALL( readInstance(scip, filenameX, &X, &absmaxX, printinfo) );
   SCIP_CALL( readInstance(scip, filenameY, &Y, NULL, printinfo) );
   SCIP_CALL( SCIPcreateModel(scip, X, Y, &ub, lb, absmaxX, FALSE) );

   // solve the problem ...
   SCIP_CALL( SCIPsolveRCproblem(scip, X, Y, ub, lb, absmaxX, relaxfilename, timeLimit) );

   // free data
   SCIP_CALL( SCIPfreeModel(scip, X, Y) );

   SCIP_CALL( SCIPfreeTransform(scip) );
   SCIP_CALL( SCIPfree(&scip) );

   BMScheckEmptyMemory();

   //SCIPprintMemoryDiagnostic(scip);
   return SCIP_OKAY;
}


/** main function for solving the partitioning problem */
int main(int argc, const char** argv)
{
   // check parameters
   std::vector<std::string> mainArgs;
   std::vector<std::string> otherArgs;

   std::vector<std::string> knownOptions{"s", "t", "m", "n", "d"};

   if ( ! readOptions(argc, argv, 2, 5, mainArgs, otherArgs, knownOptions) )
   {
      std::cerr << "usage: " << argv[0] << " <X file> <Y file> <ub> <lb> <name> [-s <settings>] [-t <time limit>] [-m <memory limit>] [-n <node limit>] ";
      std::cerr << "[-d <disp. freq>]" << std::endl;
      exit(1);
   }

   // get optional arguments
   const char* settings = 0;
   double timeLimit = 1e20;
   double memLimit = 1e20;
   SCIP_Longint nodeLimit = SCIP_LONGINT_MAX;
   int dispFreq = INT_MAX;
   for (long unsigned int j = 0; j < otherArgs.size(); ++j)
   {
      if ( otherArgs[j] == "s" )
         settings = otherArgs[++j].c_str();
      else if ( otherArgs[j] == "t" )
         timeLimit = atof(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "m" )
         memLimit = atof(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "n" )
         nodeLimit = atol(otherArgs[++j].c_str());
      else if ( otherArgs[j] == "d" )
         dispFreq = atoi(otherArgs[++j].c_str());
      else
         dispFreq = 1000;
   }

   // check for Y-file
   if ( mainArgs[1] == "" )
   {
      std::cerr << "Need to specify <Y file> containing the points to be separated from X." << std::endl;
      exit(1);
   }

   int ub = INT_MAX;
   if ( mainArgs.size() >= 3 && mainArgs[2] != "" )
      ub = atoi(mainArgs[2].c_str());

   int lb = 0;
   if ( mainArgs.size() >= 4 && mainArgs[2] == "" )
   {
      std::cerr << "Need to specify lower bound <lb> on RC." << std::endl;
      exit(1);
   }
   else if ( mainArgs.size() >= 4 )
      lb = atoi(mainArgs[3].c_str());

   if ( lb > ub )
   {
      std::cerr << "Lower bound " << lb << " is greater than upper bound " << ub << "." << std::endl;
      exit(1);
   }

   std::string relaxfilename = "";
   char* filename = 0;
   if ( mainArgs.size() == 5 )
   {
      char tmpname[SCIP_MAXSTRLEN];

      relaxfilename = mainArgs[4].c_str();
      strcpy(tmpname, relaxfilename.c_str());

      filename = tmpname;
   }

   // run clustering code
   SCIP_RETCODE retcode;
   retcode = computeRC(mainArgs[0], mainArgs[1], ub, lb, filename, settings, timeLimit, memLimit, nodeLimit, dispFreq);

   if ( retcode != SCIP_OKAY )
   {
      SCIPprintError(retcode);
      return -1;
   }
   return 0;
}
