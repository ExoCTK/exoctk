#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "_xsects.h"

static int find_index(const double *arr, const double elem, int low, int high);

/*Function: _interpolate_cross_sections
----------------------------------------
Input:
--XSects_out, RXsects_out: uninitialized array on which cross-sections will be 
--computed, Xsects_out represents the absorption cross section for each gas at the
--particular values of Tavg and Pavg. It is 3D and of size (numlevels-1, numbins,
--numgases). RXsects_out is the Rayleigh scattering array of size 
--(numlevels-1, numbins)

--Xsects_in is 4D array of size (len(Tgrid), len(Pgrid), len(wnogrid), numgases)
--It contains log10 of the values of cross-sections computed at each given (T,P,wno, gas)

--Tgrid, Pgrid, wnogrid: grids on which XSects_in was computed

--Pavg, Tavg: 1D grids of size (numlevels -1) on which to compute 
cross-sections. 

--wno_array: 1D grid of size numbins on which to compute cross-sections

--scatter_power: Rayleigh scattering power
--scatter_coeff: Rayleigh scattering coefficient
--numlevels, numbins, numsources, Pgridlen, Tgridlen: integer size of arrays
--binoffset: starting index in wnogrid for which we wish to calculate cross-sections

This function interpolates the values of the cross-sections at 
each level, wavenumber and gas by using the multivariate linear interpolation formula
(such as seen in Wikipedia: Multivariate Interpolation).  It interpolates
the log10 of the cross-sections vs the log10 of the pressure and temperature at
each level, wavenumber and gas.
*/

void _init_xsects(double *Xsects_out, double *RXsects_out, const double *Xsects_in, 
		  const double *Tgrid, const double *Pgrid, const double *Tavg, const double *Pavg,
		  const double *wno_array, double scatter_coeff, double scatter_power, 
		  int numlevels, int numbins, int binoffset, int nv, int numgases, int Tgridlen,
		  int Pgridlen) {
    
  int i, j, k;
  double sigma0 = 2.3E-27*scatter_coeff;
  double wno0 = 1E4/0.43;
  int totalbins = nv;
  int wnomin = binoffset;
  int wnomax = binoffset + numbins;
  int offset1 = Pgridlen*Tgridlen*totalbins;

  //apply log10 to each element of Pavg, Tavg, Pgrid, Tgrid so
  //we can interpolate on a log-log grid
  double *log10Pavg = malloc((numlevels-1) * sizeof(double));
  double *log10Tavg = malloc((numlevels-1) * sizeof(double));
  double *log10Tgrid = malloc(Tgridlen * sizeof(double));
  double *log10Pgrid = malloc(Pgridlen * sizeof(double));
  for (i=0; i < numlevels-1; i++) {
    log10Pavg[i] = log10(Pavg[i]);
    log10Tavg[i] = log10(Tavg[i]);
  }
  for (i=0; i < Pgridlen; i++) {
    log10Pgrid[i] = log10(Pgrid[i]);
  }
  for (i=0; i < Tgridlen; i++) {
    log10Tgrid[i] = log10(Tgrid[i]);
  }

  for (k=0; k < numlevels-1; k++) {
    double T = log10Tavg[k];
    double P = log10Pavg[k];
    int temploc = find_index(log10Tgrid, T, 0, Tgridlen - 1);
    int presloc = find_index(log10Pgrid, P, 0, Pgridlen - 1);

    double T1 = log10Tgrid[temploc];
    double T2 = log10Tgrid[temploc + 1];
    double P1 = log10Pgrid[presloc];
    double P2 = log10Pgrid[presloc+1];
    
    double totalweight = 1/((T2-T1)*(P2-P1));
    
    /*
    calculate the offsets once in this loop rather than
    numbins*numgases times in following the for loops
    while 4 dimensional arrays do exist in c, the way numpy stores
    arrays is like a large 1D array so we must index the arrays
    as if they are a 1D array
    */

    int pres_offset1 = presloc*Tgridlen*totalbins;
    int pres_offset2 = (presloc+1)*Tgridlen*totalbins;
    int temp_offset1 = temploc*totalbins;
    int temp_offset2 = (temploc+1)*totalbins;
    int level_offset_abs = k*numbins*numgases;
    double weight11 = (T2-T)*(P2-P);
    double weight21 = (T-T1)*(P2-P);
    double weight12 = (T2-T)*(P-P1);
    double weight22 = (T-T1)*(P-P1);
    for (i=wnomin; i < wnomax; i++) {
      int wno_offset = i;
      for (j=0; j < numgases; j++) {
	int gas_offset = j*offset1;
	double sigma11 = Xsects_in[gas_offset + pres_offset1 + temp_offset1 + wno_offset];
	double sigma12 = Xsects_in[gas_offset + pres_offset2 + temp_offset1 + wno_offset];
	double sigma21 = Xsects_in[gas_offset + pres_offset1 + temp_offset2 + wno_offset];
	double sigma22 = Xsects_in[gas_offset + pres_offset2 + temp_offset2 + wno_offset];
	double sigmainterp = totalweight * (sigma11*weight11 + sigma21*weight21
					    + sigma12*weight12 + sigma22*weight22);
	//Remember to exponentiate sigmainterp since we interpolated log10(sigma)
	Xsects_out[level_offset_abs + (i-binoffset)*numgases + j] = pow(10,sigmainterp);
      }
    }
    //Create Rayleigh scattering cross-sections, no interpolation here
    int level_offset_scat = k*numbins;
    for (i=0; i < numbins; i++) {
      RXsects_out[level_offset_scat + i] = sigma0*pow(wno_array[i]/wno0, 
					      scatter_power)*1E-4;
    }
  }

  //Remember to free up the memory we malloc'd at the beginning
  free(log10Tavg);
  free(log10Pavg);
  free(log10Tgrid);
  free(log10Pgrid);
}

/*
Function: find_index
--------------------
searches array for the closest value of elem. If the element is bounded by the 
array such that min(arr) < elem < max(arr), then this function returns the index i such that
arr[i] <= elem < arr[i+1]. Otherwise if elem > max(arr), it returns len(arr) - 2 and if 
elem < min(arr), it returns 0.
Uses a recursive binary search to find the index.
*/
static int find_index(const double *arr, const double elem, int low, int high) {
  //Base Case: this is where the element would be found
  if (low == high || low == (high - 1)) return low;
  //Otherwise look at the half of the array that the element
  //should be found
  int mid = (high + low) / 2;
  if (elem < arr[mid]) return find_index(arr, elem, low, mid);
  else return find_index(arr, elem, mid, high);
}

