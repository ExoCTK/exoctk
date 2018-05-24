#include <math.h>
#include "_tau.h"
#include <stdio.h>

/*Function: _calc_tau
--------------------
Input:
--Xsects: array of absorption cross sections in form [numlevels-1, numbins, numsources]
--RXsects: array of rayleigh scattering cross sections if form
  numlevels - 1 x numbins
--Z: array of vertical heights of atm layers, size: numlevels
--Pavg: array of heights of layers, size: numlevels - 1
--Tavg: array of temperatures of layers, size: numlevels - 1
--Fractions: array of fractional abundance of each species, assumed to be independent of layer level
--numlevels, numbins, numsoucres: size of arrays
--Tau: will be modified by _tau to hold each tau as a function of layer and a function of wavenumber. size: numbins x numlevels
--kb, r0: constants
*/

void _calc_tau(const double *Xsects, const double *RXsects, const double *Z, const double *Pavg, const double* Tavg, const double *Fractions, int numlevels, int numbins, int numsources, double *Tau, double kb, double r0) {
  
  int i,j,k,g;
  //printf("%d\n",numlevels);
  for (i = 0; i < numlevels - 2; i++) {
    //printf("%d %f\n",i,Fractions[i]);
    for (j = 0; j < i; j++) {
      double r1 = r0 + Z[i];
      double r2 = r0 + Z[i-j];
      double r3 = r0 + Z[i-j-1];
      double dx = sqrt(r3*r3 - r1*r1) - sqrt(r2*r2 - r1*r1);
      int curlevel = i-j-1;
      
      int offset = curlevel*numbins*numsources;
      double mult = dx*1E5*Pavg[curlevel]
	/ (kb*Tavg[curlevel]);

	//Rest of the species
      //for (g = 2; g < numsources; g ++) {
      //  printf("%d %d %.10e\n",g, curlevel, Fractions[g*(numlevels)+curlevel]);
        //printf("%d\n", g);
      //}
      for (k=0; k< numbins; k++) {
	//printf("%f\n",Fractions[curlevel]);
	double xsec;
	double dtau = 0;
	int index_offset = offset + k*numsources;
	//Rayleigh Scattering
	xsec = RXsects[curlevel*numbins + k];
	dtau += Fractions[0+curlevel]*xsec;

	//H2
	xsec = Xsects[index_offset + 0];
	dtau += Fractions[0+curlevel]*Fractions[0+curlevel]*xsec;
	//printf("%f\n",Fractions[curlevel]);
	//He
	xsec = Xsects[index_offset + 1];
	dtau += Fractions[0+curlevel]*Fractions[1*numlevels+curlevel]*xsec;
	
	//Rest of the species
	for (g = 2; g < numsources; g ++) {
	  xsec = Xsects[index_offset + g];
	  dtau += Fractions[g*numlevels+curlevel]*xsec;
          //printf("%d\n", g);
	}
	Tau[k*numlevels + i] += 2*dtau*mult;
      }
    }
  }
}
