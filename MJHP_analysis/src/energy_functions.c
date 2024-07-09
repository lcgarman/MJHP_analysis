#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"

void find_energy_bounds(Wavefunction * WFK, EnergyStep * ESTP) 
{
  int k;
  int n;
  int nkpt;
  int nband;
  double fermi;
  double bandE_min_Ha;
  double bandE_max_Ha;
  int bandE_min_integer;
  int bandE_max_integer;
  double zero;

  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi = WFK->fermi;

  printf( "\nFinding Energy Range to Scan.\n");
  printf( "\tnkpt=%d\tnband=%d\tfermi=%lf\n", nkpt, nband, fermi);
  /*find the minimum and maximum band energy*/
  bandE_min_Ha = 0.0;
  bandE_max_Ha = 0.0;
  for (k=0;k<nkpt;k++) {
    for (n=0;n<nband;n++) {
      if (WFK->eigen[k][n] < bandE_min_Ha) bandE_min_Ha = WFK->eigen[k][n];
      if (WFK->eigen[k][n] > bandE_max_Ha) bandE_max_Ha = WFK->eigen[k][n];
    }
  }
  
  /*convert E Ha min/max to eV and store*/
  ESTP->bandE_min = (bandE_min_Ha-fermi)*HATOEV;  
  ESTP->bandE_max = (bandE_max_Ha-fermi)*HATOEV;
  if (ESTP->bandE_max > 10.0) {
    ESTP->bandE_max = 10.0;
  }
  bandE_min_integer = floor(ESTP->bandE_min);
  bandE_max_integer = ceil(ESTP->bandE_max);
  
  /*find dE that corresponds to 0 eV or at the Ef*/
  zero = -(ESTP->bandE_min*EMESH)+0.5;
  ESTP->dE_zero = floor(zero);

  /*find the number of step from Elow to Ehigh using Energy mesh of 100*/
  ESTP->nEstep = (EMESH*bandE_max_integer)-(EMESH*bandE_min_integer)+1;
  printf( "\tScanning Energy Range: %d eV to %d eV (%d steps)\n", bandE_min_integer, bandE_max_integer, ESTP->nEstep);
  printf( "\tdE Step Corresponding to 0 eV (E_F) = %d\n", ESTP->dE_zero);
 
} //END of find_energy_bounds function
  

  
