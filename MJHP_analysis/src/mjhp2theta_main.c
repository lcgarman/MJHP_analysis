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
#include "energy_functions.h"
#include "psp_functions.h"
#include "twotheta_functions.h"
#include "mjhp2theta_local_potential.h"
#include "mjhp2theta_nonlocal_potential.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  char MJINfilename [200];
  char MJOUTfilename [200];
  char RFLCfilename [200];
  char ABOfilename [200];
  char ABOUTfilename [200];
  char POTfilename [200];
  char DENfilename [200];
  char WFKfilename [200];
  double fermi;
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
  MottJonesConditions mjc;
  TwoTheta tth;
  EnergyContribution econ;

  if (argc > 1) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
  }
  else {
    printf("USAGE: mjhp2theta_analyze <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_Pseudopotential(&psp);
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);
  Initialize_EnergyContribution(&econ);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);
  /*store HKL indices as conventional */

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(MJOUTfilename, fcab.MJOUTfilename);
  printf( "mjin: %s\n", ABOfilename);
  printf( "mjout: %s\n", MJOUTfilename);
  
  /*read reflection file*/
  strcpy(RFLCfilename, MJOUTfilename);
  strcat(RFLCfilename, ".rflc");
  printf( "rflc: %s\n", RFLCfilename);
  printf( "\nReading Reflections from %s\n", RFLCfilename);
  read_reflections(RFLCfilename, &tth); 

  /*Read _i_Density for Fermi energy file*/
  strcpy(DENfilename, ABOfilename);
  strcat(DENfilename, "_i_DEN");
  read_binary_abinit(DENfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom);
  fermi = wave.fermi;

  /*Read Potential file*/
  strcpy(POTfilename, ABOfilename);
  strcat(POTfilename, "_o_POT");
  read_binary_abinit(POTfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom);

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(&ucell, &grid);
  
  /*Perform FFT on Real space potential grid*/
  FFTon_RealGrid(&bin, &grid, &ucell);

  /*Read Wavefunction file*/
  strcpy(WFKfilename, ABOfilename);
  strcat(WFKfilename, "_o_WFK");
  read_binary_abinit(WFKfilename, 0, &ucell, &grid, &symm, &wave, &bin, &atom); 
  printf("\n");

  /*Find band energy range to scan*/
  wave.fermi = fermi;
  find_energy_bounds(&wave, &estp); 

  /*Find nonsymmorphic Symmetry related to HKL points*/
  symmorphic_symmetry(&symm); 
  
  /*calculate the local potential energy contribution*/
  mjhp2theta_local_potential(&tth, &estp, &grid, &wave, &ucell, &bin, &econ);

  /*Read psp information from *out file*/
  strcpy(ABOUTfilename, ABOfilename);
  strcat(ABOUTfilename, ".out");
  read_PSPdata(ABOUTfilename, &psp, &atom); 

  /*calculate the nonlocal potential energy contribution*/
  mjhp2theta_nonlocal_potential(&tth, &estp, &wave, &psp, &ucell, &atom, &econ);

  /*combine nonlocal and local potential energy*/
  concatinate_twotheta_potential(&econ, &estp, &tth); 
  
  /*calculate the Fermi sphere angle*/
  Calculate_FermiDegree(&ucell, &fsph);
 
  /*print total potential energy contributions*/ 
  strcat(MJOUTfilename, "_2theta");
  strcat(MJOUTfilename, ".mjout");
  printf( "\nPrinting potential energy contributions to %s\n", MJOUTfilename);
  print_mjhp2theta_energy(MJOUTfilename, &tth, &estp, &ucell, &fcab, &econ, &fsph);

  /*Free allocated memory not needed anymore*/
  printf( "\nFree Allocated Variables\n");
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

  tth.rflc_mult = FreeMemory_oneD_int(tth.rflc_mult);
  tth.two_theta = FreeMemory_oneD_double(tth.two_theta);
  tth.BZkpt = FreeMemory_twoD_double(tth.BZkpt, tth.nrflc);
  tth.hpw = FreeMemory_oneD_int(tth.hpw);
  tth.kpw = FreeMemory_oneD_int(tth.kpw);
  tth.lpw = FreeMemory_oneD_int(tth.lpw);
  tth.rflc_H = FreeMemory_oneD_int(tth.rflc_H);
  tth.rflc_K = FreeMemory_oneD_int(tth.rflc_K);
  tth.rflc_L = FreeMemory_oneD_int(tth.rflc_L);

  printf( "\n\nEND OF FILE\n");

  return(0);
}
