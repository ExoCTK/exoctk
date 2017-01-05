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

/*------------------- main_transmission.c ------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

/* Main procedure for computing transmission spectra  

------------------------------------------------------------------ */

#include <stdio.h>

#include "opac.h"
#include "atmos.h"
#include "prototypes.h"


/* --- Global variables ------------------------------------------ */

struct Atmos atmos;
struct Opac opac;
struct Chem chem;

/* ------- begin ---------------- main --------------------------- */


void transmission(struct vars variables, double* wavelength, double* flux)
{
  TotalOpac(variables);
  printf("TotalOpac done\n");

  ReadTP(variables);
  printf("ReadTP done\n");

  RT_Transmit(variables, wavelength, flux);
  printf("Exo_Transmit done\n");

  FreeTP(variables);
  FreeOpacTable(variables, opac);
  FreeChemTable(variables);
}

/* ------- end -------- main_transmission.c  --------------------- */
