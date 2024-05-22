#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "rwa_functions.h"

void read_mjin_header(char filename[100], MottJonesConditions * MJC, FileCabinet* FCAB, FermiSphere* FS)
{
  FILE* mjin_file;
  char strv[20];
  char ABOfilename[100];
  char MJOUTfilename[100];
  char lattice[10];
  double vec; 
  int jzH, jzK, jzL;
  double scanE_low, scanE_high;

  printf("Reading %s\n", filename);
  mjin_file=fopen(filename,"r+"); /*open the filename in read/write mode*/
  if(mjin_file==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  fscanf(mjin_file, "%s\n", ABOfilename);
  strcpy(FCAB->ABOfilename, ABOfilename);
  fscanf(mjin_file, "%s\n", MJOUTfilename);
  strcpy(FCAB->MJOUTfilename, MJOUTfilename);
  fscanf(mjin_file, "%s", strv);
  fscanf(mjin_file, "%lf", &vec);
  FS->vec = vec;
  fscanf(mjin_file, "%s", lattice);
  strcpy(MJC->lattice, lattice);

  fscanf(mjin_file, "%s", strv);
  fscanf(mjin_file, "%d", &jzH);
  fscanf(mjin_file, "%d", &jzK);
  fscanf(mjin_file, "%d", &jzL);
  MJC->jzH = jzH;
  MJC->jzK = jzK;
  MJC->jzL = jzL;

  fscanf(mjin_file, "%s", strv);
  fscanf(mjin_file, "%lf", &scanE_low);
  fscanf(mjin_file, "%lf", &scanE_high);
  MJC->scanE_low = scanE_low;
  MJC->scanE_high = scanE_high;

  fclose(mjin_file);
} //END of readSS 

void print_mjhpHKL_energy(char filename[200], VectorIndices *VECT, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON)
{
  FILE* fhkl;
  int nEstep;
  int dE;
  double Elow;
  int H, K, L;
  double total_potential;

  nEstep = ESTP->nEstep;
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;

  fhkl = fopen(filename, "w");
  if(fhkl==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  printf("\nPrinting Energy for HKL to %s\n", filename);

  /*print header information*/
  fprintf(fhkl, "%s\n", FCAB->ABOfilename);
  fprintf(fhkl, "astar: %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);  
  fprintf(fhkl, "bstar: %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);  
  fprintf(fhkl, "cstar: %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);  
  fprintf(fhkl, "noSteps %d \n", ESTP->nEstep);
  Elow = ESTP->bandE_min*EMESH;
  fprintf(fhkl, "lowestE %lf\n", floor(Elow));
  fprintf(fhkl, "Emesh %d\n", EMESH);
  fprintf(fhkl, "HKL %d\t%d\t%d\n", H, K, L);
  fprintf(fhkl, "\n");
  /*END of Header Info*/

  /*Print Energy by hkl*/
  total_potential = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    fprintf(fhkl, "%d\t%e\n", dE, ECON->total_potential[dE]);
    total_potential += ECON->total_potential[dE];
  }
  printf("Total Potential Energy for %d %d %d = %e\n", VECT->H, VECT->K, VECT->L, total_potential);

  /*free energy contribution*/
  ECON->total_potential = FreeMemory_oneD_double(ECON->total_potential);
} //END of Print_Reflection function

void print_mjhp2theta_energy(char filename[200], TwoTheta *TTH, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON, FermiSphere *FS)
{
  FILE* f2th;
  int nEstep;
  int dE;
  int nrflc;
  int n;
  double Elow;
  double reflection_total;

  nEstep = ESTP->nEstep;
  nrflc = TTH->nrflc;

  f2th = fopen(filename, "w");
  if(f2th==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  printf("\nPrinting Energy for MJHP_twotheta to %s\n", filename);

  /*print header information*/
  fprintf(f2th, "%s\n", FCAB->ABOfilename);
  fprintf(f2th, "astar: %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);  
  fprintf(f2th, "bstar: %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);  
  fprintf(f2th, "cstar: %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);  
  fprintf(f2th, "noSteps %d \n", ESTP->nEstep);
  Elow = ESTP->bandE_min*EMESH;
  fprintf(f2th, "lowestE %lf\n", floor(Elow));
  fprintf(f2th, "Emesh %d\n", EMESH);
  fprintf(f2th, "FS_angle %lf\n", FS->two_theta);
  fprintf(f2th, "nrflc %d\n", nrflc);
  fprintf(f2th, "\n");
  /*END of Header Info*/

  /*Print Energy by hkl*/
  for (n=0;n<nrflc;n++) {
    fprintf(f2th, "\n");
    fprintf(f2th, "rflc#%d\t %lf\t %d %d %d\t%d\n", n, TTH->two_theta[n], TTH->rflc_H[n], TTH->rflc_K[n], TTH->rflc_L[n], TTH->rflc_mult[n]);
    reflection_total = 0.0;
    for (dE=0;dE<nEstep;dE++) {
      fprintf(f2th, "%d\t%e\n", dE, ECON->rflc_total[dE][n]);
      reflection_total += ECON->rflc_total[dE][n];
    }
  }

  /*free energy contribution*/
  ECON->rflc_total = FreeMemory_twoD_double(ECON->rflc_total, nEstep);
} //END of Print_Reflection function

