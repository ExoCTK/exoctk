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

/*--------------------------- lin_interpol.c ---------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

#include <stdio.h>
#include <math.h>
#include <rpc/types.h>

#include "prototypes.h"

/* ------- begin ---------------- lint2D.c ----------------------- */

/* 2-D Linear interpolation procedure:  Finds the point on a 
   plane (x,y,z) that lies between these grid points: 

   * (x1,y1,z1)   * (x2,y1,z2)

   * (x1,y2,z3)   * (x2,y2,z4)
------------------------------------------------------------------ */

double lint2D (double x1, double x2, double y1, double y2, double z1,  double z2, double z3, double z4, double x, double y)
{
  double z_int1, z_int2, z;
  
  z_int1 = lint(x1, z1, x2, z2, x);
  z_int2 = lint(x1, z3, x2, z4, x);
  z = lint(y1, z_int1, y2, z_int2, y);

  return z;
}

/* ------- end ------------------ lint2D.c ----------------------- */

/* ------- begin ---------------- lint.c ------------------------- */

/* 1-D Linear interpolation procedure:  Finds the point on a 
   line (x,y) that lies between (x1,y1) and (x2,y2).
------------------------------------------------------------------ */

double lint(double x1, double y1, double x2, double y2, double x)
{
  double slope, intcept, y;
  slope = (y2 - y1) / (x2 - x1);
  intcept = y1 - slope * x1;
  y = intcept + slope * x;
  return y;
}

/* ------- end ------------------ lint.c ------------------------- */


