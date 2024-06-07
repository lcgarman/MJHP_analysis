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
  FILE* flog;
  char MJINfilename [200];
  char MJOUTfilename [200];
  char MJLOGfilename [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char POTfilename [200];
  char WFKfilename [200];
  char H_append [10];
  char K_append [10];
  char L_append [10];
  int H_conventional, K_conventional, L_conventional;
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
  /*store HKL indices as conventional */
  H_conventional = mjc.jzH;
  K_conventional = mjc.jzK;
  L_conventional = mjc.jzL;

  /*open log file and print header*/
  strcpy(MJLOGfilename, fcab.MJOUTfilename);
  strcat(MJLOGfilename, "_");
  sprintf(H_append, "%d", H_conventional);
  strcat(MJLOGfilename, H_append);
  sprintf(K_append, "%d", K_conventional);
  strcat(MJLOGfilename, K_append);
  sprintf(L_append, "%d", L_conventional);
  strcat(MJLOGfilename, L_append);
  strcat(MJLOGfilename, ".mjlog");
  flog = fopen(MJLOGfilename, "w");
  if(flog==NULL) {
    printf("%s not found. \n", MJLOGfilename);
    exit(0);
  }  

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJOUTfilename, fcab.MJOUTfilename);
  fprintf(flog, "mjin: %s\n", ABOfilename);
  fprintf(flog, "mjout: %s\n", MJOUTfilename);
  /*transfrom HKL to primitive centering*/
  transform_HKL(&mjc); 
  /*store primitive HKL in VectorIndices array*/
  vect.H = mjc.jzH;
  vect.K = mjc.jzK;
  vect.L = mjc.jzL;
  printf("H=%d K=%d L=%d (%d %d %d)\n", vect.H, vect.K, vect.L, H_conventional, K_conventional, L_conventional); 
  fprintf(flog, "\nHKL conventional: %d %d %d\n", H_conventional, K_conventional, L_conventional);
  fprintf(flog, "HKL primitive: %d %d %d\n", vect.H, vect.K, vect.L); 

  /*Read Potential file*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  read_binary_abinit(flog, POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom);

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(flog, &ucell, &grid);
  
  /*Perform FFT on Real space potential grid*/
  FFTon_RealGrid(flog, &bin, &grid, &ucell);

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  read_binary_abinit(flog, WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom); 
  printf("\n");

  /*Find band energy range to scan*/
  find_energy_bounds(flog, &wave, &estp); 

  /*Find nonsymmorphic Symmetry related to HKL points*/
  symmorphic_symmetry(&symm); 
  find_symmetric_hkl(flog, &vect, &symm, &grid); 
  
  /*find the shell in reciprocal space to include potential energy contributions*/
  find_MJregion(flog, &vect, &ucell);
  
  /*calculate the local potential energy contribution*/
  mjhpHKL_local_potential(flog, &vect, &grid, &estp, &wave, &ucell, &bin, &econ);

  /*Read psp information from *out file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  read_PSPdata(flog, ABOUTfilename, &psp, &atom); 

  /*calculate the nonlocal potential energy contribution*/
  mjhpHKL_nonlocal_potential(flog, &vect, &grid, &estp, &wave, &psp, &ucell, &atom, &econ);

  /*combine nonlocal and local potential energy*/
  concatinate_HKL_potential(flog, &econ, &estp);

  /*Integrate the total potential energy up to Ef*/
  integrate_HKL_potential(flog, &econ, &estp, &atom);
  
  /*print total potential energy contributions*/ 
  strcat(MJOUTfilename, "_");
  strcat(MJOUTfilename, H_append);
  strcat(MJOUTfilename, K_append);
  strcat(MJOUTfilename, L_append);
  strcat(MJOUTfilename, ".mjout");
  fprintf(flog, "\nPrinting potential energy contributions to %s\n", MJOUTfilename);
  print_mjhpHKL_energy(MJOUTfilename, &vect, &estp, &ucell, &fcab, &econ);

  /*Free allocated memory not needed anymore*/
  fprintf(flog, "\nFree Allocated Variables\n");
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
  symm.symor = FreeMemory_threeD_int(symm.symor, 3, 3);
  vect.H_arr = FreeMemory_oneD_int(vect.H_arr);
  vect.K_arr = FreeMemory_oneD_int(vect.K_arr);
  vect.L_arr = FreeMemory_oneD_int(vect.L_arr);
  vect.H_posarr = FreeMemory_oneD_int(vect.H_posarr);
  vect.K_posarr =  FreeMemory_oneD_int(vect.K_posarr);
  vect.L_posarr = FreeMemory_oneD_int(vect.L_posarr);

  fprintf(flog, "\n\nEND OF FILE\n");
  fclose(flog);

  return(0);
}
