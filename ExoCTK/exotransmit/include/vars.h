/* This file is part of Exo_transmit.

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

/*----------------------- getVars.h ---------------------------------

Author: Patrick Slough (sloughpa@grinnell.edu)

--------------------------------------------------------------------- */

#ifndef getVars_h
#define getVars_h

/* --- Vars structure ----------------------------------------------- */

typedef struct vars{
  int NTAU;
  
  int NTEMP;
  double TLOW;
  double THIGH;
  
  int NPRESSURE;
  double PLOW;
  double PHIGH;
  double THRESHOLD;
  double RAYLEIGH;
  
  int NLAMBDA;
  double G;
  double R_PLANET;
  double R_STAR;
  double T_STAR;

  char* tpfname;
  char* eosfname;

  int chemselection[32];

}vars;

//extern vars variables;

/* --- Associated function prototypes ------------------------------- */

vars getVars();
void getNTau(vars *variables, char *filename);

#endif

/* ------- end ------------------------- getVars.h ------------------ */
