/* This file is part of Exo_transmit.

    Exo_transmit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Exo_transmit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Exo_Transmit.  If not, see <http://www.gnu.org/licenses/>.
*/

/*-------------------------- geometry.c --------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

/* Routines dealing with the geometry for calculating a 
   transmission spectrum. 

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "opac.h"
#include "atmos.h"
#include "constant.h"
#include "include.h"
#include "prototypes.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Opac opac;

/* ------- begin ---------------- Radius ------------------------- */

/* Given the planet radius, this determines the radius
for the planet plus atmosphere: R.

------------------------------------------------------------------ */

double Radius(double R_pl, double *ds, int NT)
{
  double R;
  int i;
  
  R = R_pl;

  for(i=0; i<NT; i++){
    R += ds[i];
  }

  return R;
}

/* ------- end ------------------ Radius ------------------------- */

/* ------- begin ---------------- Angles ------------------------- */

/* Determines theta and dtheta values for the flux integral.

------------------------------------------------------------------ */

void Angles(double *ds, double *theta, double *dtheta, int NT, double Rpl, double Rst)
{
  double h, R;
  int j;

  R = Radius(Rpl, ds, NT);
  
  h = R;
  printf("R (planet & atmosphere) %e\n", h);

  for(j=0; j<NT; j++){

    h -= ds[j];
    theta[j] = asin(h/Rst);

    if(j==0)
      dtheta[j] = asin(R/Rst) - theta[j];
    else
      dtheta[j] = theta[j-1] - theta[j];

  }

  return;  
}

/* ------- end ------------------ Angles ------------------------- */

/* ------- begin ---------------- Tau_LOS ------------------------ */

/* Determines line-of-sight optical depths tau_tr as a function 
of wavelength and height in the atmosphere.

------------------------------------------------------------------ */

void Tau_LOS(double **kappa_nu, double **tau_tr, double *ds, int NT, double Rpl, int NLam)
{
  double R, a, b, dl[NT];
  int i, j, k;

  R = Radius(Rpl, ds, NT);
  printf("R_atmosphere %e\n", R - Rpl);

  for(j=0; j<NT; j++)
    dl[j] = 0.0;

  a = R;

  for(j=0; j<NT; j++){

    a -= ds[j];
    b = R;

    for(k=0; k<=j; k++){
      dl[k] = 2.0 * pow(SQ(b) - SQ(a), 0.5);
      b -= ds[k];
    }

    for(k=0; k<=j; k++){
      if (k != j)
	dl[k] -= dl[k+1];
    }
    
    for(k=0; k<=j; k++){
      for(i=0; i<NLam; i++){
	
	/* Calculates tau_tr by summing up kappa * dl 
	   for shells along the line of sight */
	tau_tr[i][j] += kappa_nu[i][k] * dl[k];

      }
    }
  }
  
  return;
}

/* ------- end ------------------ Tau_LOS ------------------------ */
