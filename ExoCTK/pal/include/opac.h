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

/*--------------------------- opac.h --------------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

--------------------------------------------------------------------- */

#ifndef __OPAC_H__
#define __OPAC_H__

/* --- Opacity structure -------------------------------------------- */

struct Opac { 
	char* name;
  int NP, NT;
  double ***kappa, **abundance;
  double *P, *Plog10, *T;
};

/* --- Chemistry structure ------------------------------------------ */

struct Chem {
  double **total, **C, **CH4, **CO, **CO2, **C2H2, **C2H4, **C2H6, **H, 
    **HCN, **HCl, **HF, **H2, **H2CO, **H2O, **H2S, **He, **K, **MgH, 
    **N, **N2, **NO2, **NH3, **NO, **Na, **O, **O2, **O3, **OH, **OCS, 
    **PH3, **SH, **SO2, **SiH, **SiO, **TiO, **VO, **mu;
  double *P, *T;
};

/* --- Associated function prototypes ------------------------------- */

void ReadOpacTable(struct Opac opac, char *filename);
void FreeOpacTable(struct Opac opac);

#endif /* !__OPAC_H__ */

/* ------- end ---------------------------- opac.h ------------------ */
