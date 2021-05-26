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
