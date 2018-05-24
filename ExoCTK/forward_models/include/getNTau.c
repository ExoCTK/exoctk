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

/*-------------------------- getNTau.c ---------------------------

Author: Patrick Slough (sloughpa@grinnell.edu)

------------------------------------------------------------------ */

/* Determines how many layers of the T-P profile lie above the 
pressure of the cloud layer specified in userInput.in

------------------------------------------------------------------ */


#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vars.h"

void getNTau(vars *variables, char *filename){
  if(variables->THRESHOLD != 0){
    char dum[7]; 
    int NTau = 0;
    int i = 0;
    double pressure[variables->NTAU];
    int index[variables->NTAU];
    double threshold = variables->THRESHOLD;
    
    printf("Cloud layer at %lf pa\n", threshold);
    
    FILE *file = fopen(filename, "r");
    if(file){
      fscanf(file, "%s %s %s", dum, dum, dum);
      while(fscanf(file, "%d %lf %*lf", &index[i], &pressure[i]) != EOF){
	if(pressure[i] < threshold)
	  NTau++;
	i++;
      }
      NTau++;
    }
    else{
      printf("Error reading File\n");
      exit(1);
    }
    
    fclose(file);
    
    variables->NTAU = NTau;
  }
  else
    printf("No Cloud layer will be calculated\n");
}

/* -------- end ----------------- GetNTau ------------------------ */

