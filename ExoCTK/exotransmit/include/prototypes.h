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

/*----------------------- protytopes.h ------------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

--------------------------------------------------------------------- */


#ifndef __PROTOTYPES_H__
#define __PROTOTYPES_H__

/*---- Function prototypes associated w/ Exo_Transmit --------------- */

void getChemSelection(int array[]);
char **getFileArray();
double Radius(double R_pl, double *ds, int NT);
void Angles(double *ds, double *theta, double *dtheta, int NT, double Rpl, double Rst);
void Tau_LOS(double **kappa_nu, double **tau_tr, double *ds, int NT, double Rpl, int NLam);
void Locate(int n, double *array, double value, int *ilow);
int  RT_Transmit(struct vars variables, double* wavelength, double* flux);
void ReadTP(vars variables);
void TotalOpac(struct vars variables);
double Planck(double T, double lambda);
double lint2D(double x1, double x2, double y1, double y2, double z1, 
	      double z2, double z3, double z4, double x, double y);
double lint(double x1, double y1, double x2, double y2, double x);
void FreeTP(struct vars variables);
void ReadChemTable(struct vars variables);
void FreeChemTable(struct vars variables);
void errorCheck(int onOff, double value);

#endif

/*---- end ---------------- prototypes.h ---------------------------- */
