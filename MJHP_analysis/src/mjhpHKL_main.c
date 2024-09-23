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
#include "fftw_functions.h"
#include "hkl_functions.h"
#include "energy_functions.h"
#include "mjhpHKL_local_potential.h"
#include "mjhpHKL_nonlocal_potential.h"
#include "psp_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char MJOUTfilename_N [200];
  char MJOUTfilename [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char POTfilename [200];
  char WFKfilename [200];
  char H_append [10];
  char K_append [10];
  char L_append [10];
  int H_conventional, K_conventional, L_conventional;
  int n;
  int ndts;
  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  EnergyStep estp;
  FermiSphere fsph; 
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  Pseudopotential psp;
  VectorIndices vect;
  MottJonesConditions mjc;
  EnergyContribution econ;

  if (argc > 1) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: mjhpHKL_analyze <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_Pseudopotential(&psp);
  Initialize_VectorIndices(&vect);
  Initialize_MottJonesConditions(&mjc);
  Initialize_EnergyContribution(&econ);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);
  ndts = mjc.ndts;

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJOUTfilename, fcab.MJOUTfilename);
  printf("mjin: %s\n", ABOfilename);
  printf("mjout: %s\n", MJOUTfilename);
  printf("\n Number of HKL: %d\n", ndts);
  for(n=0;n<ndts;n++) {
    printf("DTSET %d\t HKL = %d %d %d\n", n+1, mjc.jzH[n], mjc.jzK[n], mjc.jzL[n]);
  }
  printf("\nBegin Reading Files...\n");

  /*Read Potential file*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  read_binary_abinit(POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom);

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);
  
  /*Perform FFT on Real space potential grid*/
  FFTon_RealGrid(&bin, &grid, &ucell);
  bin.cc_rec_grid= FreeMemory_threeD_complex(bin.cc_rec_grid, grid.ngfftx, grid.ngffty);

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  read_binary_abinit(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom); 

  /*Find band energy range to scan*/
  find_energy_bounds(&wave, &estp); 

  /*Read psp information from *out file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  read_PSPdata(ABOUTfilename, &psp, &atom); 
  /*end of reading in info*/

  /*being calcing potential for HKL sets*/
  for (n=0;n<ndts;n++) {
    /*store HKL indices as conventional */
    mjc.nH = mjc.jzH[n];
    mjc.nK = mjc.jzK[n];
    mjc.nL = mjc.jzL[n];
  	H_conventional = mjc.nH; 
  	K_conventional = mjc.nK; 
  	L_conventional = mjc.nL; 
    printf( "\nDTSET %d\t HKL = %d %d %d\n", n+1, mjc.nH, mjc.nK, mjc.nL);

	/*transfrom HKL to primitive centering*/
	transform_HKL(&mjc); 
	/*store primitive HKL in VectorIndices array*/
	vect.H = mjc.nH;
	vect.K = mjc.nK;
	vect.L = mjc.nL;
	printf("H=%d K=%d L=%d (%d %d %d)\n", vect.H, vect.K, vect.L, H_conventional, K_conventional, L_conventional); 
	printf( "\nHKL conventional: %d %d %d\n", H_conventional, K_conventional, L_conventional);
	printf( "HKL primitive: %d %d %d\n", vect.H, vect.K, vect.L); 
	
	/*Find nonsymmorphic Symmetry related to HKL points*/
	find_symmetric_hkl(&vect, &symm, &grid); 
	
	/*find the shell in reciprocal space to include potential energy contributions*/
	find_MJregion(&vect, &ucell);
	
    //wave.fermi = 0.07601;
	/*calculate the local potential energy contribution*/
	mjhpHKL_local_potential(&vect, &grid, &estp, &wave, &ucell, &bin, &econ);
	
    /*calculate the nonlocal potential energy contribution*/
	mjhpHKL_nonlocal_potential(&vect, &grid, &estp, &wave, &psp, &ucell, &atom, &econ);

	/*combine nonlocal and local potential energy*/
	concatinate_HKL_potential(&econ, &estp);

	/*Integrate the total potential energy up to Ef*/
	integrate_HKL_potential(&econ, &estp, &atom);
	
	/*print total potential energy contributions*/ 
	sprintf(H_append, "%d", H_conventional);
	sprintf(K_append, "%d", K_conventional);
	sprintf(L_append, "%d", L_conventional);
    strcpy(MJOUTfilename_N, MJOUTfilename);
	strcat(MJOUTfilename_N, "_");
	strcat(MJOUTfilename_N, H_append);
	strcat(MJOUTfilename_N, K_append);
	strcat(MJOUTfilename_N, L_append);
	strcat(MJOUTfilename_N, ".mjout");
	printf( "\nPrinting potential energy contributions to %s\n", MJOUTfilename_N);
	print_mjhpHKL_energy(MJOUTfilename_N, &vect, &estp, &ucell, &fcab, &econ);
 
	/*free HKL variables for next iteration*/
	vect.H_arr = FreeMemory_oneD_int(vect.H_arr);
	vect.K_arr = FreeMemory_oneD_int(vect.K_arr);
	vect.L_arr = FreeMemory_oneD_int(vect.L_arr);
	vect.H_posarr = FreeMemory_oneD_int(vect.H_posarr);
	vect.K_posarr =  FreeMemory_oneD_int(vect.K_posarr);
	vect.L_posarr = FreeMemory_oneD_int(vect.L_posarr);
  }


  /*Free allocated memory not needed anymore*/
  printf( "\nFree Allocated Variables\n");
  FreeMemory_Wavefunctions(&wave);
  bin.rec_grid = FreeMemory_threeD_complex(bin.rec_grid, grid.ngfftx, grid.ngffty);
  FreeMemory_PSPvariables(&psp);
  wave.npw = FreeMemory_oneD_int(wave.npw);
  symm.symrel = FreeMemory_threeD_int(symm.symrel, 3, 3);
  atom.typat = FreeMemory_oneD_int(atom.typat);
  wave.kpt = FreeMemory_twoD_double(wave.kpt, 3);
  symm.tnons = FreeMemory_twoD_double(symm.tnons, 3);
  wave.wtk = FreeMemory_oneD_double(wave.wtk);
  symm.mult = FreeMemory_oneD_int(symm.mult);
  atom.xred = FreeMemory_twoD_double(atom.xred, 3);
  wave.eigen = FreeMemory_twoD_double(wave.eigen, wave.nkpt);
  mjc.jzH = FreeMemory_oneD_int(mjc.jzH);
  mjc.jzK = FreeMemory_oneD_int(mjc.jzK);
  mjc.jzL = FreeMemory_oneD_int(mjc.jzL);

  printf( "\n\nEND OF FILE\n");

  return(0);
}
