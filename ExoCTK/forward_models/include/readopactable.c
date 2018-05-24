/* This file is part of Exo_Transmit.

    Exo_Transmit is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Exo_Transmit is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Exo_Transmit.  If not, see <http://www.gnu.org/licenses/>.
*/

/*----------------------- readopactable.c ------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include "atmos.h"
#include "opac.h"
#include "nrutil.h" 
#include "constant.h" 
#include "vars.h"
#include "prototypes.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
 
/* ---------------------------------------------------------------
 * Read in opacity files: kappa(pressure, temperature, lambda)
 * Opacity file actually contains cross section in units of m^2.  
 * Need to multiply by the number density to get kappa in units of
 * m^-1
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadOpacTable.c -------------------- */

void ReadOpacTable(struct Opac opac, char *filename) {

  int i, j, k;
  double junk;
  vars variables = getVars();
  int NLAMBDA = variables.NLAMBDA;
  int NTEMP = variables.NTEMP;
  int NPRESSURE = variables.NPRESSURE;
  
  FILE *f1;

  atmos.lambda = dvector(0, NLAMBDA-1);

  opac.NP = NPRESSURE;
  opac.NT = NTEMP; 
 
  f1 = fopen(filename,"r");
  if(f1 == NULL){
    printf("\nreadopactable.c:\nError opening %s opacity file\nPlease check that the proper name and path is specified in otherInput.in\n", opac.name);
    exit(1);
  }
  
  for (k=0; k<opac.NT; k++) {
    fscanf(f1,"%le", &opac.T[k]);
  }
  
  for (j=0; j<opac.NP; j++) {
    fscanf(f1,"%le", &opac.P[j]);
    opac.Plog10[j] = log10(opac.P[j]);
  }
  
  for (i=0; i<NLAMBDA; i++) {
    fscanf(f1,"%le", &atmos.lambda[i]);
	
    for (j=0; j<opac.NP; j++) {
      fscanf(f1,"%le", &junk);
	  	  
      for (k=0; k<opac.NT; k++) {
		  fscanf(f1,"%le", &opac.kappa[i][j][k]);
	
	/*   kappa in file is actually cross section, sigma.  
	     Need to multiply by number density */
	
	opac.kappa[i][j][k] *= opac.abundance[j][k] * opac.P[j] /
	  (KBOLTZMANN * opac.T[k]);
	
      }
    }
  }
  
  fclose(f1);
  printf("opac %e %e %e\n", atmos.lambda[NLAMBDA-1], opac.P[0], opac.T[0]);
    
}


/* ------- end -------------- ReadOpacTable.c -------------------- */

/* ------- begin ------------ FreeOpacTable.c -------------------- */

void FreeOpacTable(struct Opac opac)
{
  vars variables = getVars();
  int NLAMBDA = variables.NLAMBDA;
  int NTEMP = variables.NTEMP;
  int NPRESSURE = variables.NPRESSURE;

  free_dvector(opac.T, 0, NTEMP-1);
  free_dvector(opac.P, 0, NPRESSURE-1);
  free_dvector(opac.Plog10, 0, NPRESSURE-1);
  free_dmatrix(opac.abundance, 0, NPRESSURE-1, 0, NTEMP-1);
  free_d3tensor(opac.kappa, 0, NLAMBDA-1, 0, NPRESSURE-1, 0, NTEMP-1);

}

/* ------- end -------------- FreeOpacTable.c -------------------- */
