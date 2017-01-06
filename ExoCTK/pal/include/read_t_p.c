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


/*------------ file ------- read_t_p.c ---------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

/* Reads in the temperature - pressure profile from the file 
   defined in userIinput.in
------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "atmos.h"
#include "nrutil.h"
#include "opac.h"
#include "prototypes.h"
#include "vars.h"

/* --- Global variables ------------------------------------------ */

extern struct Atmos atmos;
extern struct Chem chem;
extern struct Opac opac;

/* ------- begin --------------------- ReadTP.c ------------------ */

void ReadTP(struct vars variables)
{
  
  /* Get relevant variables */
  getNTau(&variables, variables.tpfname);
  
  /* Rename for convenience */
  int NTEMP = variables.NTEMP; 
  int NPRESSURE = variables.NPRESSURE;
  int NTAU = variables.NTAU;
  int THRESHOLD = variables.THRESHOLD;
  
  int i, j, k, num[NTAU];
  char dum[7];
  FILE *file;
  
  /* Allocate memory for atmos structure */
  
  atmos.P = dvector(0, NTAU-1);
  atmos.T = dvector(0, NTAU-1);
  atmos.mu = dvector(0, NTAU-1);
  
  file = fopen(variables.tpfname, "r");						
  if(file == NULL){
    printf("\nread_t_p.c:\nError opening file: No such file or directory\n\n");
    exit(1);
  }
  
  /* Read in T-P profile */

  fscanf(file, "%s %s %s", dum, dum, dum);
  for (i=0; i<NTAU; i++){
    fscanf(file, "%d %le %le", &num[i],  &atmos.P[i], &atmos.T[i]);
  }
  fclose(file);
  
  /* Determine mean molecular weight at each altitude */
  
  /* Interpolation variables */
  int x1,x2,y1,y2;
  double t1,t2,t,p1,p2,p,z1,z2,z3,z4;
  for (i=0;i<NTAU;i++){ 			
    
    /* locate temperature on grid */    
    j = 0;
    while((atmos.T[i] > chem.T[j]) && (atmos.T[i] > chem.T[j+1]) 
	  && (j <= (NTEMP-1)))     
      j++;
    if (j > (NTEMP-1))
      printf("Tempearature not found on grid\n");
    else{
      /* set interpolation points*/
      x1 = j;
      x2 = j+1;
      t1 = chem.T[j]; 
      t2 = chem.T[j+1];
      t = atmos.T[i];
    }
    
    /* locate pressure on grid */
    k=0;
    while((atmos.P[i] > chem.P[k]) && (atmos.P[i] >= chem.P[k+1]) 
	  && (k <= (NPRESSURE-1)))
      k++;
    if (k > (NPRESSURE-1))
      printf("Pressure not found on grid\n");
    else{
      /* set interpolation points*/
      y1 =k;
      y2 = k+1;
      p1 = chem.P[y1]; 
      p2 = chem.P[y2];
      p = atmos.P[i];
    }
    
    /* weighted atmospheric mass at (p1,t1) (p1,t2) (p2,t1) (p2,t2) */
    z1 = chem.mu[y1][x1];
    z2 = chem.mu[y1][x2];
    z3 = chem.mu[y2][x1];
    z4 = chem.mu[y2][x2];
    
    /* linear interpolate */
    atmos.mu[i] = lint2D(t1,t2,p1,p2,z1,z2,z3,z4,t,p);
  }
  
  /* Cloud layer calculation */

  if(THRESHOLD != 0.0){			 
    double proportion = (log10(THRESHOLD) - log10(atmos.P[NTAU-2]))
      / (log10(atmos.P[NTAU-1]) - log10(atmos.P[NTAU-2]));  
    atmos.P[NTAU-1] = THRESHOLD;
    atmos.T[NTAU-1] = proportion*(atmos.T[NTAU-1] - atmos.T[NTAU-2]) 
      + atmos.T[NTAU-2];
    atmos.mu[NTAU-1] = proportion*(atmos.mu[NTAU-1] - atmos.mu[NTAU-2]) 
      + atmos.mu[NTAU-2];
  }
  
  
  printf("%d\t%e\t%e\n", num[NTAU-1], atmos.P[NTAU-1], atmos.T[NTAU-1]);
    
}

/* ------- end ----------------------- ReadTP.c ------------------ */

/* ------- start --------------------- FreeTP.c ------------------ */

void FreeTP(struct vars variables){
  
  /* Frees atmos structure */
  
  // vars variables = getVars();
  int NTAU = variables.NTAU;
  int NLAMBDA = variables.NLAMBDA;

  free_dvector(atmos.P, 0, NTAU-1);
  free_dvector(atmos.T, 0, NTAU-1);
  free_dvector(atmos.mu, 0, NTAU-1);
  free_dvector(atmos.lambda, 0, NLAMBDA-1);

}

/* ------- end ----------------------- FreeTP.c ------------------ */
