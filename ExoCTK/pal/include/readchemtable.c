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

/*----------------------- readchemtable.c ------------------------

Author: Eliza Kempton (kemptone@grinnell.edu)

------------------------------------------------------------------ */

#include <stdio.h>
#include <stdlib.h>

#include "opac.h" 
#include "nrutil.h" 
#include "vars.h"
#include "prototypes.h"
 
/* --- Global variables ------------------------------------------ */

extern struct Chem chem;

/* ---------------------------------------------------------------
 * Read in chemistry files: abundance(pressure, temperature)
 * --------------------------------------------------------------- */

/* ------- begin ------------ ReadChemTable.c -------------------- */

void ReadChemTable(struct vars variables) {
  
  int i, j;
  char dum[8];
  int chemSelection[32];
  
  /* Initialize and obtain variables from other files */
  // char **fileArray = getFileArray();
  // vars variables = getVars();
  // getChemSelection(chemSelection);
  for (i=0; i<32; i++) {
    chemSelection[i] = variables.chemselection[i];
  }
  

  int NTEMP = variables.NTEMP;
  int NPRESSURE = variables.NPRESSURE;
  
  printf("NTEMP: %d, NPRESSURE: %d\n", NTEMP, NPRESSURE);		
  
  FILE *f1;
  
  /* Allocate memory for Chem structure */
  
  chem.T = dvector(0, NTEMP-1);	
  chem.P = dvector(0, NPRESSURE-1);
  
  chem.total = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.C = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.CH4 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.CO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.CO2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.C2H2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.C2H4 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.C2H6 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);		
  chem.H = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.HCN = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.HCl = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.HF = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.H2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.H2CO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.H2O = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.H2S = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.He = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.K = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.MgH = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.N = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.N2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.NO2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.NH3 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.NO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.Na = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);
  chem.O = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.O2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.O3 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.OCS = dmatrix(0, NPRESSURE-1, 0, NTEMP-1); 	
  chem.OH = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.PH3 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.SH = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.SO2 = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.SiH = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.SiO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.TiO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  chem.VO = dmatrix(0, NPRESSURE-1, 0, NTEMP-1);	
  
  /* Read in chemistry table */	
  
  // f1 = fopen(fileArray[1],"r");
  printf("Using EOS file: %s\n", variables.eosfname);
  f1 = fopen(variables.eosfname, "r");
  if(f1 == NULL){
    // printf("\nreadchemtable.c:\nError opening file: %s \nNo such file or directory.\nMake sure you have the appropriate file path and name specified in userInput.in\n\n", fileArray[1]);
    exit(1);
  }
  
  for (i=0; i<38; i++)
    fscanf(f1,"%s", dum);
  
  for (i=NPRESSURE-1; i>=0; i--){
    fscanf(f1,"%le", &chem.P[i]);
    for (j=NTEMP-1; j>=0; j--) {
      if(i == NPRESSURE-1){
	fscanf(f1,"%le", &chem.T[j]); 
	fscanf(f1,"%le", &chem.total[i][j]);
	fscanf(f1,"%le", &chem.C[i][j]);
	fscanf(f1,"%le", &chem.CH4[i][j]); 
	/* Checks for errors, but only the first time */
	errorCheck(chemSelection[0], chem.CH4[i][j]);
	fscanf(f1,"%le", &chem.CO[i][j]);  
	errorCheck(chemSelection[2], chem.CO[i][j]);
	fscanf(f1,"%le", &chem.OCS[i][j]); 
	errorCheck(chemSelection[19], chem.OCS[i][j]);
	fscanf(f1,"%le", &chem.CO2[i][j]); 
	errorCheck(chemSelection[1], chem.CO2[i][j]);
	fscanf(f1,"%le", &chem.C2H2[i][j]); 
	errorCheck(chemSelection[7], chem.C2H2[i][j]);
	fscanf(f1,"%le", &chem.C2H4[i][j]); 
	errorCheck(chemSelection[8], chem.C2H4[i][j]);
	fscanf(f1,"%le", &chem.C2H6[i][j]); 
	errorCheck(chemSelection[8], chem.C2H6[i][j]);
	fscanf(f1,"%le", &chem.H[i][j]);
	fscanf(f1,"%le", &chem.HCN[i][j]); 
	errorCheck(chemSelection[13], chem.HCN[i][j]);
	fscanf(f1,"%le", &chem.HCl[i][j]); 
	errorCheck(chemSelection[12], chem.HCl[i][j]);
	fscanf(f1,"%le", &chem.HF[i][j]); 
	errorCheck(chemSelection[14], chem.HF[i][j]);
	fscanf(f1,"%le", &chem.H2[i][j]);
	fscanf(f1,"%le", &chem.H2CO[i][j]); 
	errorCheck(chemSelection[10], chem.H2CO[i][j]);
	fscanf(f1,"%le", &chem.H2O[i][j]); 
	errorCheck(chemSelection[3], chem.H2O[i][j]);
	fscanf(f1,"%le", &chem.H2S[i][j]); 
	errorCheck(chemSelection[11], chem.H2S[i][j]);
	fscanf(f1,"%le", &chem.He[i][j]);
	fscanf(f1,"%le", &chem.K[i][j]);
	fscanf(f1,"%le", &chem.MgH[i][j]); 
	errorCheck(chemSelection[15], chem.MgH[i][j]);
	fscanf(f1,"%le", &chem.N[i][j]);
	fscanf(f1,"%le", &chem.N2[i][j]); 
	errorCheck(chemSelection[16], chem.N2[i][j]);
	fscanf(f1,"%le", &chem.NO2[i][j]); 
	errorCheck(chemSelection[18], chem.NO2[i][j]); 
	fscanf(f1,"%le", &chem.NH3[i][j]); 
	errorCheck(chemSelection[4], chem.NH3[i][j]);
	fscanf(f1,"%le", &chem.NO[i][j]); 
	errorCheck(chemSelection[17], chem.NO[i][j]);
	fscanf(f1,"%le", &chem.Na[i][j]); 
	fscanf(f1,"%le", &chem.O[i][j]);
	fscanf(f1,"%le", &chem.O2[i][j]); 
	errorCheck(chemSelection[5], chem.O2[i][j]);
	fscanf(f1,"%le", &chem.O3[i][j]); 
	errorCheck(chemSelection[5], chem.O3[i][j]);
	fscanf(f1,"%le", &chem.OH[i][j]); 
	errorCheck(chemSelection[20], chem.OH[i][j]);
	fscanf(f1,"%le", &chem.PH3[i][j]); 
	errorCheck(chemSelection[21], chem.PH3[i][j]); 
	fscanf(f1,"%le", &chem.SH[i][j]); 
	errorCheck(chemSelection[22], chem.SH[i][j]);
	fscanf(f1,"%le", &chem.SO2[i][j]); 
	errorCheck(chemSelection[25], chem.SO2[i][j]); 
	fscanf(f1,"%le", &chem.SiH[i][j]); 
	errorCheck(chemSelection[23], chem.SiH[i][j]); 
	fscanf(f1,"%le", &chem.SiO[i][j]); 
	errorCheck(chemSelection[24], chem.SiO[i][j]); 
	fscanf(f1,"%le", &chem.TiO[i][j]); 
	errorCheck(chemSelection[26], chem.TiO[i][j]); 
	fscanf(f1,"%le", &chem.VO[i][j]); 
	errorCheck(chemSelection[27], chem.VO[i][j]);
      }
      else{
	fscanf(f1,"%le", &chem.T[j]); 
	fscanf(f1,"%le", &chem.total[i][j]);
	fscanf(f1,"%le", &chem.C[i][j]);
	fscanf(f1,"%le", &chem.CH4[i][j]); 
	fscanf(f1,"%le", &chem.CO[i][j]); 
	fscanf(f1,"%le", &chem.OCS[i][j]); 
	fscanf(f1,"%le", &chem.CO2[i][j]); 
	fscanf(f1,"%le", &chem.C2H2[i][j]); 
	fscanf(f1,"%le", &chem.C2H4[i][j]); 
	fscanf(f1,"%le", &chem.C2H6[i][j]); 
	fscanf(f1,"%le", &chem.H[i][j]);
	fscanf(f1,"%le", &chem.HCN[i][j]); 
	fscanf(f1,"%le", &chem.HCl[i][j]); 
	fscanf(f1,"%le", &chem.HF[i][j]); 
	fscanf(f1,"%le", &chem.H2[i][j]);
	fscanf(f1,"%le", &chem.H2CO[i][j]); 
	fscanf(f1,"%le", &chem.H2O[i][j]); 
	fscanf(f1,"%le", &chem.H2S[i][j]); 
	fscanf(f1,"%le", &chem.He[i][j]);
	fscanf(f1,"%le", &chem.K[i][j]);
	fscanf(f1,"%le", &chem.MgH[i][j]); 
	fscanf(f1,"%le", &chem.N[i][j]);
	fscanf(f1,"%le", &chem.N2[i][j]); 
	fscanf(f1,"%le", &chem.NO2[i][j]); 
	fscanf(f1,"%le", &chem.NH3[i][j]); 
	fscanf(f1,"%le", &chem.NO[i][j]); 
	fscanf(f1,"%le", &chem.Na[i][j]); 
	fscanf(f1,"%le", &chem.O[i][j]);
	fscanf(f1,"%le", &chem.O2[i][j]); 
	fscanf(f1,"%le", &chem.O3[i][j]); 
	fscanf(f1,"%le", &chem.OH[i][j]); 
	fscanf(f1,"%le", &chem.PH3[i][j]); 
	fscanf(f1,"%le", &chem.SH[i][j]); 
	fscanf(f1,"%le", &chem.SO2[i][j]); 
	fscanf(f1,"%le", &chem.SiH[i][j]); 
	fscanf(f1,"%le", &chem.SiO[i][j]); 
	fscanf(f1,"%le", &chem.TiO[i][j]); 
	fscanf(f1,"%le", &chem.VO[i][j]); 
      }
    }
  }
  
  fclose(f1);
  
  /* Print out final line of data as a double-check */
  printf("Chemistry: \n");	
  printf("P_0\t%e \n", chem.P[NPRESSURE-1]);
  printf("T_0\t%e \n", chem.T[NTEMP-1]);
  printf("total\t%e \n", chem.total[0][0]);
  printf("C\t%e \n", chem.C[0][0]);
  printf("CH4\t%e \n", chem.CH4[0][0]);
  printf("CO\t%e \n", chem.CO[0][0]);
  printf("CO2\t%e \n", chem.CO2[0][0]);
  printf("C2H2\t%e \n", chem.C2H2[0][0]);
  printf("C2H4\t%e \n", chem.C2H4[0][0]);
  printf("H\t%e \n", chem.H[0][0]);
  printf("HCN\t%e \n", chem.HCN[0][0]);
  printf("HCl\t%e \n", chem.HCl[0][0]);
  printf("HF\t%e \n", chem.HF[0][0]);
  printf("H2\t%e \n", chem.H2[0][0]);
  printf("H2CO\t%e \n", chem.H2CO[0][0]);
  printf("H2O\t%e \n", chem.H2O[0][0]);
  printf("H2S\t%e \n", chem.H2S[0][0]);
  printf("He\t%e \n", chem.He[0][0]);
  printf("MgH\t%e \n", chem.MgH[0][0]);
  printf("N\t%e \n", chem.N[0][0]);
  printf("N2\t%e \n", chem.N2[0][0]);
  printf("NO2\t%e \n", chem.NO2[0][0]);
  printf("NH3\t%e \n", chem.NH3[0][0]);
  printf("NO\t%e \n", chem.NO[0][0]);
  printf("O\t%e \n", chem.O[0][0]);
  printf("O2\t%e \n", chem.O2[0][0]);
  printf("OCS\t%e \n", chem.OCS[0][0]);
  printf("OH\t%e \n", chem.OH[0][0]);
  printf("PH3\t%e \n", chem.PH3[0][0]);
  printf("SH\t%e \n", chem.SH[0][0]);
  printf("SO2\t%e \n", chem.SO2[0][0]);
  printf("SiH\t%e \n", chem.SiH[0][0]);
  printf("SiO\t%e \n", chem.SiO[0][0]);
  printf("TiO\t%e \n", chem.TiO[0][0]);
  printf("VO\t%e \n", chem.VO[0][0]);
  printf("O3\t%e \n", chem.O3[0][0]);
  printf("C2H6\t%e \n", chem.C2H6[0][0]);
  printf("Na\t%e \n", chem.Na[0][0]);
  printf("K\t%e \n", chem.K[0][0]);
  
  return;
}

/* ------- end -------------- ReadChemTable.c -------------------- */

/* ------- start ------------ FreeChemTable.c -------------------- */

void FreeChemTable(struct vars variables){ 

  // vars variables = getVars();
  int NTEMP = variables.NTEMP;
  int NPRESSURE = variables.NPRESSURE;

  free_dvector(chem.T, 0, NTEMP-1);	
  free_dvector(chem.P, 0, NPRESSURE-1);
  
  free_dmatrix(chem.total, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.C, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.CH4, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.CO, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.CO2, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.C2H2, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.C2H4, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.C2H6, 0, NPRESSURE-1, 0, NTEMP-1);		
  free_dmatrix(chem.H, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.HCN, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.HCl, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.HF, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.H2, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.H2CO, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.H2O, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.H2S, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.He, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.K, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.MgH, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.N, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.N2, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.NO2, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.NH3, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.NO, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.Na, 0, NPRESSURE-1, 0, NTEMP-1);
  free_dmatrix(chem.O, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.O2, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.O3, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.OCS, 0, NPRESSURE-1, 0, NTEMP-1); 	
  free_dmatrix(chem.OH, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.PH3, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.SH, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.SO2, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.SiH, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.SiO, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.TiO, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.VO, 0, NPRESSURE-1, 0, NTEMP-1);	
  free_dmatrix(chem.mu, 0, NPRESSURE-1, 0, NTEMP-1);	

}

/* ------- end -------------- FreeChemTable.c -------------------- */

/* ------- start ------------ errorCheck.c ----------------------- */

/* Shows error message if gas not selected is present in chem table. */
void errorCheck(int onOff, double value){
  if(onOff == 0 && value != 0.0){
    printf("Error: Gas not selected for inclusion in Calculation has non-zero Chem table entry\n");
  }
}

/* ------- end -------------- errorCheck.c ----------------------- */
