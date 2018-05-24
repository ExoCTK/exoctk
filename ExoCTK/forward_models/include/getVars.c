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

/*-------------------------- getVars.c ---------------------------

Author: Patrick Slough (sloughpa@grinnell.edu)

------------------------------------------------------------------ */

/* Gets input variables from userInput.in

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "vars.h"

vars getVars(){
  
  vars vars;
  FILE *userVarInput = fopen("userInput.in", "r");
  
  int lineCount = 1;
  
  if(userVarInput){
    char line[256];
    while(fgets(line, sizeof(line), userVarInput) != NULL){
      switch(lineCount){
      case 12:
	vars.G = strtod(line, NULL);
	lineCount++;
	break;
      case 14:
	vars.R_PLANET = strtod(line, NULL);
	lineCount++;
	break;
      case 16:
	vars.R_STAR = strtod(line, NULL);
	lineCount++;
	break;
       case 18:
	vars.THRESHOLD = strtod(line, NULL);
	lineCount++;
	break;
      case 20:
	vars.RAYLEIGH = strtod(line, NULL);
	lineCount++;
	break;
      default:
	lineCount++;
	break;
      }
    }
  }
  else{
    printf("Error opening file\n");
    exit(1);
  }

  fclose(userVarInput);
  
  FILE *otherVarInput = fopen("otherInput.in", "r");
  
  lineCount = 1;
  
  if(otherVarInput){
    char line[256];
    while(fgets(line, sizeof(line), otherVarInput) != NULL){
      switch(lineCount){
      case 38:
	vars.NTEMP = strtol(line, NULL, 10);
	lineCount++;
	break;
      case 39:
	vars.TLOW = strtod(line, NULL);
	lineCount++;
	break;
      case 40:
	vars.THIGH = strtod(line, NULL);
	lineCount++;
	break;
      case 42:
	vars.NPRESSURE = strtol(line, NULL, 10);
	lineCount++;
	break;
      case 43:
	vars.PLOW = strtod(line, NULL);
	lineCount++;
	break;
      case 44:
	vars.PHIGH = strtod(line, NULL);
	lineCount++;
	break;
      case 46:
	vars.NLAMBDA = strtol(line, NULL, 10);
	lineCount++;
	break;
      case 48:
	vars.NTAU = strtol(line, NULL, 10);
	lineCount++;
	break;
      default:
	lineCount++;
	break;
      }
    }
  }
  else{
    printf("Error opening file\n");
    exit(1);
  }

  fclose(otherVarInput);
  
  return vars;
  
}         


/* -------- end ----------------- GetVars ------------------------ */

