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


/*-------------------------- getChemSelection.c ------------------

Author: Patrick Slough (sloughpa@grinnell.edu)

------------------------------------------------------------------ */

/* Determines which molecules will be included in the opacity 
   calculation based on user input from selectChem.in

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ------- begin ---------------- GetChemSelection --------------- */

void getChemSelection(int array[]){
  
  FILE *input;
  int i, lineCount;

  i = 0;
  lineCount = 1;

  input = fopen("selectChem.in", "r");
  
  if(input){
    char line[256];
    while(fgets(line, sizeof(line), input) != NULL){
      if(lineCount >= 2 && lineCount < 34){
	int index = 0;
	while(line[index] != '\n'){
	  if(line[index] == '1' || line[index] == '0')
	    break;
	  else
	    index++;
	}
	array[i++] = atoi(&line[index]);
	lineCount++;
      }
      else{
	lineCount++;
      }
    }
  }
  else{
    printf("Error reading input file\n");
    exit(1);
  }

  fclose(input);

}

/* -------- end ----------------- GetChemSelection --------------- */
