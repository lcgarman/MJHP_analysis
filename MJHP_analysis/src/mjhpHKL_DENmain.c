#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "read_binary_abinit.h"
#include "cell_parameters.h"
#include "rwa_functions.h"
#include "hkl_functions.h"
#include "mjhpHKL_density.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char MJXSFfilename [200];
  char MJXSFfilename_N [200];
  char MJXSFfilename_N1 [200];
  char MJXSFfilename_N2 [200];
  char ABOfilename [200];
  char WFKfilename [200];
  char H_append [10];
  char K_append [10];
  char L_append [10];
  char scanEmin_append [20];
  char scanEmid_append [20];
  char scanEmax_append [20];
  int H_conventional, K_conventional, L_conventional;
  int n;
  int ndts;
  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  FermiSphere fsph;
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  VectorIndices vect;
  MottJonesConditions mjc;

  if (argc > 1) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: mjhpHKL_density <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_VectorIndices(&vect);
  Initialize_MottJonesConditions(&mjc);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);
  ndts = mjc.ndts;
  sprintf(scanEmin_append, "%lf", mjc.scanE_min);
  sprintf(scanEmid_append, "%lf", mjc.scanE_mid);
  sprintf(scanEmax_append, "%lf", mjc.scanE_max);

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJXSFfilename, fcab.MJOUTfilename);
  printf("mjin: %s\n", ABOfilename);
  printf("mjXSF: %s\n", MJXSFfilename);
  printf("\n Number of HKL: %d\n", ndts);
  for(n=0;n<ndts;n++) {
    printf("DTSET %d\t HKL = %d %d %d\n", n+1, mjc.jzH[n], mjc.jzK[n], mjc.jzL[n]);
  }
  printf("\nBegin Reading Files...\n");

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  read_binary_abinit(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom); 

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);

  for (n=0;n<ndts;n++) {
    mjc.nH = mjc.jzH[n];
    mjc.nK = mjc.jzK[n];
    mjc.nL = mjc.jzL[n];
	/*store HKL indices as conventional */
	H_conventional = mjc.nH;
	K_conventional = mjc.nK;
	L_conventional = mjc.nL;
    printf("\nDTSET %d\t HKL = %d %d %d\n", n+1, mjc.nH, mjc.nK, mjc.nL);
    printf("\n\nDTSET %d\t HKL = %d %d %d\n", n+1, mjc.nH, mjc.nK, mjc.nL);

	/*transfrom HKL to primitive centering*/
	transform_HKL(&mjc); 
	/*store primitive HKL in VectorIndices array*/
	vect.H = mjc.nH;
	vect.K = mjc.nK;
	vect.L = mjc.nL;
	printf("H=%d K=%d L=%d (%d %d %d)\n", vect.H, vect.K, vect.L, H_conventional, K_conventional, L_conventional); 
	printf("\nHKL conventional: %d %d %d\n", H_conventional, K_conventional, L_conventional);
	printf("HKL primitive: %d %d %d\n", vect.H, vect.K, vect.L); 
	
    /*find symmetric HKL indices*/
	find_symmetric_hkl(&vect, &symm, &grid); 
	
	/*find the shell in reciprocal space to include potential energy contributions*/
	find_MJregion(&vect, &ucell);
	
	/*prepare xsffilename_HKL_*/
	sprintf(H_append, "%d", H_conventional);
	sprintf(K_append, "%d", K_conventional);
	sprintf(L_append, "%d", L_conventional);
    strcpy(MJXSFfilename_N, MJXSFfilename);
	strcat(MJXSFfilename_N, "_");
	strcat(MJXSFfilename_N, H_append);
	strcat(MJXSFfilename_N, K_append);
	strcat(MJXSFfilename_N, L_append);
	strcat(MJXSFfilename_N, "_");
	strcpy(MJXSFfilename_N1, MJXSFfilename_N);
	strcpy(MJXSFfilename_N2, MJXSFfilename_N);
	strcat(MJXSFfilename_N1, scanEmin_append);
	strcat(MJXSFfilename_N1, "_");
	strcat(MJXSFfilename_N1, scanEmid_append);
	strcat(MJXSFfilename_N1, ".xsf");
	strcat(MJXSFfilename_N2, scanEmid_append);
	strcat(MJXSFfilename_N2, "_");
	strcat(MJXSFfilename_N2, scanEmax_append);
	strcat(MJXSFfilename_N2, ".xsf");
	
	/* FIRST calc density from minE to 0*/
	vect.scanE_start = mjc.scanE_min;
	vect.scanE_stop = mjc.scanE_mid;
	printf("\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
	printf("\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
	/*calculate the local potential energy contribution*/
	mjhpHKL_density(&grid, &bin, &wave, &symm, &ucell, &vect);
	/*print total potential energy contributions*/ 
	printf("\nPrinting Mott-Jones Density to: %s\n", MJXSFfilename_N1);
	print_XSF(MJXSFfilename_N1, &ucell, &grid, &bin, &atom);
	
	/* SECOND calc density from 0 maxE*/
	vect.scanE_start = mjc.scanE_mid;
	vect.scanE_stop = mjc.scanE_max;
	printf("\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
	printf("\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
	/*calculate the local potential energy contribution*/
	mjhpHKL_density(&grid, &bin, &wave, &symm, &ucell, &vect);
	/*print total potential energy contributions*/ 
	printf("\nPrinting Mott-Jones Density to: %s\n", MJXSFfilename_N2);
	print_XSF(MJXSFfilename_N2, &ucell, &grid, &bin, &atom);

	/*free HKL variables for next iteration*/
	vect.H_arr = FreeMemory_oneD_int(vect.H_arr);
	vect.K_arr = FreeMemory_oneD_int(vect.K_arr);
	vect.L_arr = FreeMemory_oneD_int(vect.L_arr);
	vect.H_posarr = FreeMemory_oneD_int(vect.H_posarr);
	vect.K_posarr =  FreeMemory_oneD_int(vect.K_posarr);
	vect.L_posarr = FreeMemory_oneD_int(vect.L_posarr);
  }
	
  /*Free allocated memory not needed anymore*/
  printf("\nFree Allocated Variables\n");
  FreeMemory_Wavefunctions(&wave);
  wave.npw = FreeMemory_oneD_int(wave.npw);
  symm.symrel = FreeMemory_threeD_int(symm.symrel, 3, 3);
  atom.typat = FreeMemory_oneD_int(atom.typat);
  wave.kpt = FreeMemory_twoD_double(wave.kpt, 3);
  symm.tnons = FreeMemory_twoD_double(symm.tnons, 3);
  wave.wtk = FreeMemory_oneD_double(wave.wtk);
  symm.mult = FreeMemory_oneD_int(symm.mult);
  atom.xred = FreeMemory_twoD_double(atom.xred, 3);
  wave.eigen = FreeMemory_twoD_double(wave.eigen, wave.nkpt);

  printf("\n\nEND OF FILE\n");

  return(0);
}
