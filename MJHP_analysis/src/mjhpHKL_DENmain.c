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
  FILE* flog;
  char MJINfilename [200];
  char MJXSFfilename [200];
  char MJXSFfilename1 [200];
  char MJXSFfilename2 [200];
  char MJLOGfilename [200];
  char ABOfilename [200];
  char WFKfilename [200];
  char H_append [10];
  char K_append [10];
  char L_append [10];
  char scanEmin_append [10];
  char scanEmid_append [10];
  char scanEmax_append [10];
  int H_conventional, K_conventional, L_conventional;
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
  /*store HKL indices as conventional */
  H_conventional = mjc.jzH;
  K_conventional = mjc.jzK;
  L_conventional = mjc.jzL;

  /*open log file and print header*/
  strcpy(MJLOGfilename, fcab.ABOfilename);
  strcat(MJLOGfilename, "_");
  sprintf(H_append, "%d", H_conventional);
  strcat(MJLOGfilename, H_append);
  sprintf(K_append, "%d", K_conventional);
  strcat(MJLOGfilename, K_append);
  sprintf(L_append, "%d", L_conventional);
  strcat(MJLOGfilename, L_append);
  strcat(MJLOGfilename, "_");
  sprintf(scanEmin_append, "%.2f", mjc.scanE_min);
  strcat(MJLOGfilename, scanEmin_append);
  strcat(MJLOGfilename, "_");
  sprintf(scanEmid_append, "%.2f", mjc.scanE_mid);
  strcat(MJLOGfilename, scanEmid_append);
  strcat(MJLOGfilename, "_");
  sprintf(scanEmax_append, "%.2f", mjc.scanE_max);
  strcat(MJLOGfilename, scanEmax_append);
  strcat(MJLOGfilename, ".mjlog");
  flog = fopen(MJLOGfilename, "w");
  if(flog==NULL) {
    printf("%s not found. \n", MJLOGfilename);
    exit(0);
  }  

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJXSFfilename, fcab.MJOUTfilename);
  fprintf(flog, "mjin: %s\n", ABOfilename);
  fprintf(flog, "mjXSF: %s\n", MJXSFfilename);
  /*transfrom HKL to primitive centering*/
  transform_HKL(&mjc); 
  /*store primitive HKL in VectorIndices array*/
  vect.H = mjc.jzH;
  vect.K = mjc.jzK;
  vect.L = mjc.jzL;
  printf("H=%d K=%d L=%d (%d %d %d)\n", vect.H, vect.K, vect.L, H_conventional, K_conventional, L_conventional); 
  fprintf(flog, "\nHKL conventional: %d %d %d\n", H_conventional, K_conventional, L_conventional);
  fprintf(flog, "HKL primitive: %d %d %d\n", vect.H, vect.K, vect.L); 

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  read_binary_abinit(flog, WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom); 
  printf("\n");

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(flog, &ucell, &grid);
  
  /*Find nonsymmorphic Symmetry related to HKL points*/
  symmorphic_symmetry(&symm); 
  find_symmetric_hkl(flog, &vect, &symm, &grid); 
  
  /*find the shell in reciprocal space to include potential energy contributions*/
  find_MJregion(flog, &vect, &ucell);
  
  /*prepare xsffilename_HKL_*/
  strcat(MJXSFfilename, "_");
  strcat(MJXSFfilename, H_append);
  strcat(MJXSFfilename, K_append);
  strcat(MJXSFfilename, L_append);
  strcat(MJXSFfilename, "_");
  strcpy(MJXSFfilename1, MJXSFfilename);
  strcpy(MJXSFfilename2, MJXSFfilename);
  strcat(MJXSFfilename1, scanEmin_append);
  strcat(MJXSFfilename1, "_");
  strcat(MJXSFfilename1, scanEmid_append);
  strcat(MJXSFfilename1, ".xsf");
  strcat(MJXSFfilename2, scanEmid_append);
  strcat(MJXSFfilename2, "_");
  strcat(MJXSFfilename2, scanEmax_append);
  strcat(MJXSFfilename2, ".xsf");

  /* FIRST calc density from minE to 0*/
  vect.scanE_start = mjc.scanE_min;
  vect.scanE_stop = mjc.scanE_mid;
  fprintf(flog, "\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
  printf("Scanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
  /*calculate the local potential energy contribution*/
  mjhpHKL_density(flog, &grid, &bin, &wave, &symm, &ucell, &vect);
  /*print total potential energy contributions*/ 
  fprintf(flog, "\nPrinting Mott-Jones Density to: %s.\n", MJXSFfilename1);
  print_XSF(MJXSFfilename1, &ucell, &grid, &bin, &atom);

  /* SECOND calc density from 0 maxE*/
  vect.scanE_start = mjc.scanE_mid;
  vect.scanE_stop = mjc.scanE_max;
  fprintf(flog, "\nScanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
  printf("Scanning Energy in Range: %lf to %lf (eV)\n", vect.scanE_start, vect.scanE_stop);
  /*calculate the local potential energy contribution*/
  mjhpHKL_density(flog, &grid, &bin, &wave, &symm, &ucell, &vect);
  /*print total potential energy contributions*/ 
  fprintf(flog, "\nPrinting Mott-Jones Density to: %s.\n", MJXSFfilename2);
  print_XSF(MJXSFfilename2, &ucell, &grid, &bin, &atom);

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
