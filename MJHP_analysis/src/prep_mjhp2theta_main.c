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
#include "rwa_functions.h"
#include "read_binary_abinit.h"
#include "cell_parameters.h"
#include "fftw_functions.h"
#include "recgrid_functions.h"
#include "twotheta_functions.h"

int main(int argc, char * argv[])
{

  /*declaring filenames*/
  FILE* flog;
  char MJINfilename [200];
  char RFLCfilename [200];
  char MJLOGfilename [200];
  char ABOfilename [200];
  char ABINfilename [200];
  char MABINfilename [200];
  char DENfilename [200];
  /*declaring structres*/
  FileCabinet fcab;
  UnitCell ucell;
  NumberGrid grid; 
  FermiSphere fsph;
  Symmetry symm;
  Wavefunction wave;
  BinaryGrid bin;
  AtomicVariables atom;
  MottJonesConditions mjc;
  TwoTheta tth;

  if (argc > 1) {
      //copy file string to inupt files
    strcpy(MJINfilename, argv[1]);
printf("working\n");
  }
  else {
    printf("USAGE: prepare_mjhp2theta <*.mjin> ");
	exit(0);
  }

  /*initialize all pointers in struct to NULL*/
  Initialize_NumberGrid(&grid);
  Initialize_Symmetry(&symm);
  Initialize_Wavefunction(&wave);
  Initialize_BinaryGrid(&bin);
  Initialize_AtomicVariables(&atom);
  Initialize_MottJonesConditions(&mjc);
  Initialize_TwoTheta(&tth);
  
  /*Read the mj file and store variables*/
  read_mjin_header(MJINfilename, &mjc, &fcab, &fsph);

  /*open log file and print header*/
  strcpy(MJLOGfilename, fcab.ABOfilename);
  strcat(MJLOGfilename, "_rflc");
  strcat(MJLOGfilename, ".mjlog");
  flog = fopen(MJLOGfilename, "w");
  if(flog==NULL) {
    printf("%s not found. \n", MJLOGfilename);
    exit(0);
  }  

  /*store file names*/
  strcpy(ABOfilename, fcab.ABOfilename);
  strcpy(RFLCfilename, fcab.MJOUTfilename);
  fprintf(flog, "mjin: %s\n", ABOfilename);
  fprintf(flog, "rflc: %s\n", RFLCfilename);

  /*read abinit *in file and modify variables*/
  strcpy(ABINfilename, ABOfilename);
  strcat(ABINfilename, ".in");
  strcpy(MABINfilename, ABINfilename);
  strcat(MABINfilename, "-modified");
  modify_abinitin(ABINfilename, MABINfilename);

  /*Read Potential file*/
  strcpy(DENfilename, ABOfilename);
  strcat(DENfilename, "_i_DEN");
  read_binary_abinit(flog, DENfilename, 1, &ucell, &grid, &symm, &wave, &bin, &atom);

  /*Find Unit Cell parameters in real and reciprocal space*/
  Determine_CellParameters(flog, &ucell, &grid);

  /*Find nonsymmorphic Symmetry */
  symmorphic_symmetry(&symm); 

  /*define the HKL grid in reciprocal space*/
  find_HKLgrid_max(&grid); 
  prepare_HKLgrid(&grid);
  
  /*Perform FFT on Real space potential grid*/
  FFTon_RealGrid(flog, &bin, &grid, &ucell);

  /*calculate the simulated powder pattern*/
  calculate_powder_pattern(flog, &tth, &bin, &grid, &ucell, &symm); 
  /*fold the reflections back into the 1st BZ*/
  fold_reflections_toBZ(flog, &tth);
  /*apply symmetry to minimize the number of kpts needed to analyze*/
  symmetry_folded_reflections(flog, MABINfilename, &tth, &symm);

  /*print total potential energy contributions*/ 
  strcat(RFLCfilename, ".rflc");
  fprintf(flog, "\nPrinting Reflection Information to %s\n", RFLCfilename);
  print_reflections(RFLCfilename, &tth);

  /*rename modified abinit in file to original in filename*/
  remove(ABINfilename);
  rename(MABINfilename, ABINfilename);

  /*Free allocated memory not needed anymore*/
  fprintf(flog, "\nFree Allocated Variables\n");
  symm.symrel = FreeMemory_threeD_int(symm.symrel, 3, 3);
  atom.typat = FreeMemory_oneD_int(atom.typat);
  wave.kpt = FreeMemory_twoD_double(wave.kpt, 3);
  symm.tnons = FreeMemory_twoD_double(symm.tnons, 3);
  symm.mult = FreeMemory_oneD_int(symm.mult);
  atom.xred = FreeMemory_twoD_double(atom.xred, 3);
  symm.symor = FreeMemory_threeD_int(symm.symor, 3, 3);

  fprintf(flog, "\n\nEND OF FILE\n");
  fclose(flog);

  return(0);
}
