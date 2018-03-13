/* ------- file: -------------------------- planck.c ----------------

       Author:        Han Uitenbroek (huitenbroek@nso.edu)
       Modified by Eliza Miller-Ricci
       Last modified: Aug 23, 2007 --

       --------------------------                      -------------- */

/* --- Evaluate the Planck function and its derivative 
       at temperature T and wavelength lambda.

       Input:  T              -- temperature [K].
               lambda         -- wavelength [m].

       Output: Bnu            -- Planck function value
                                 [J m^-2 s^-1 Hz^-1 sr^-1].
 *     --                                              -------------- */

#include <stdio.h> 
#include <math.h>

#include "constant.h"
#include "include.h"

#define MAX_EXPONENT  400.0

/* ------- begin -------------------------- Planck.c ---------------- */

double Planck(double T, double lambda)
{

  double Bnu, twohnu3_c2, hc_Tkla;

  hc_Tkla     = (HPLANCK * CLIGHT) / (T * KBOLTZMANN * lambda);
  twohnu3_c2 = (2.0*HPLANCK*CLIGHT) / CUBE(lambda); 

  if (hc_Tkla <= MAX_EXPONENT)
    Bnu = twohnu3_c2 / (exp(hc_Tkla) - 1.0); 
  else
    Bnu = 0.0;
  return Bnu;
}
/* ------- end ---------------------------- Planck.c ---------------- */


