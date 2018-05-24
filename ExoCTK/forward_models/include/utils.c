#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <rpc/types.h>

#include "include.h"
#include "nrutil.h"

/* ------------------------------------ Locate.c --------------------

       Version:       rh1.0
       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Last modified: Tue Feb 16 14:40:59 1999 --

       --------------------------                      ----------RH-- */


/* --- Find index of value in array (cf., Num. Recipes, p 90).

 Note: The Num. Recipes routine does not give the correct index
       for values that are exactly equal to an array value!
       --                                              -------------- */
 
/* ------- begin -------------------------- Locate.c ---------------- */

void Locate(int n, double *array, double value, int *ilow)
{
  bool_t ascend;
  int    ihigh, index;

  ascend = (array[n-1] > array[0]) ? TRUE : FALSE;
  *ilow = 0;  ihigh = n;

  if (ascend) {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value >= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  } else {
    while (ihigh - *ilow > 1) {
      index = (ihigh + *ilow) >> 1;
      if (value <= array[index])
	*ilow = index;
      else
	ihigh = index;
    }
  }
}
/* ------- end ---------------------------- Locate.c ---------------- */

