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
  int check;
  char lattice[10];
  double vec; 
  int n;
  int ndts;
  int jzH, jzK, jzL;
  double scanE_min, scanE_mid, scanE_max;

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
  fscanf(mjin_file, "%d", &ndts);
  MJC->ndts = ndts;
  
  /*allocate memory for HKL indices*/
  MJC->jzH = AllocateMemory_oneD_int(MJC->jzH, ndts);
  MJC->jzK = AllocateMemory_oneD_int(MJC->jzK, ndts);
  MJC->jzL = AllocateMemory_oneD_int(MJC->jzL, ndts);
 
  check=fscanf(mjin_file, "%s", strv);
  if (check==EOF) return;
  for (n=0;n<ndts;n++) {
	fscanf(mjin_file, "%d", &jzH);
	fscanf(mjin_file, "%d", &jzK);
	fscanf(mjin_file, "%d", &jzL);
	MJC->jzH[n] = jzH;
	MJC->jzK[n] = jzK;
	MJC->jzL[n] = jzL;
  }
	
  check=fscanf(mjin_file, "%s", strv);
  if (check==EOF) return;
  fscanf(mjin_file, "%lf", &scanE_min);
  fscanf(mjin_file, "%lf", &scanE_mid);
  fscanf(mjin_file, "%lf", &scanE_max);
  MJC->scanE_min = scanE_min;
  MJC->scanE_mid = scanE_mid;
  MJC->scanE_max = scanE_max;

  fclose(mjin_file);
} //END of readSS 

void print_mjhpHKL_energy(char filename[200], VectorIndices *VECT, EnergyStep *ESTP, UnitCell *UC, FileCabinet *FCAB, EnergyContribution *ECON)
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

void print_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* BIN, AtomicVariables* ATM)
{
  FILE * fxsf;
  int jx, jy, jz;
  int j;
  int natom;
  double x, y, z;
  double* Xcart;
  double* Ycart;
  double* Zcart;
  int NGX, NGY, NGZ;
  int line_counter;

  natom = ATM->natom;
  NGX = GRD->NGX;
  NGY = GRD->NGY;
  NGZ = GRD->NGZ;
  Xcart = NULL;
  Ycart = NULL;
  Zcart = NULL;

  printf("Printing MJHP Density to: %s\n", filename);
  fxsf=fopen(filename, "w");
  if (fxsf==NULL) {
    printf("ERROR: %s not found\n", filename);
    exit(0);
  }

  /*allocate mem for XYZcart*/
  Xcart = AllocateMemory_oneD_double(Xcart, natom);
  Ycart = AllocateMemory_oneD_double(Ycart, natom);
  Zcart = AllocateMemory_oneD_double(Zcart, natom);

  /*Calculate cartesian XYZ of atoms*/
  for (j=0;j<natom;j++) {
	x = ATM->xred[0][j];
	y = ATM->xred[1][j];
	z = ATM->xred[2][j];
	Xcart[j] = x*UC->ang_ax+y*UC->ang_bx+z*UC->ang_cx; 
	Ycart[j] = x*UC->ang_ay+y*UC->ang_by+z*UC->ang_cy; 
	Zcart[j] = x*UC->ang_az+y*UC->ang_bz+z*UC->ang_cz; 
  }

  fprintf(fxsf, " DIM-GROUP\n"); 
  fprintf(fxsf, " 3 1 \n");
  fprintf(fxsf, "PRIMVEC\n"); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_ax, UC->ang_ay, UC->ang_az); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_bx, UC->ang_by, UC->ang_bz); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_cx, UC->ang_cy, UC->ang_cz); 
  /* coordinates of primitive lattice */
  fprintf(fxsf, "PRIMCOORD\n" ); 
  fprintf(fxsf, "\t\t%d%3d\n", natom, 1); 
  for(j=0;j<natom;j++) {
    fprintf(fxsf, "%9d%20.10lf%20.10lf%20.10lf\n", ATM->atomicno[j], Xcart[j], Ycart[j], Zcart[j]);
  }
  fprintf(fxsf," BEGIN_BLOCK_DATAGRID3D\n");
  fprintf(fxsf," Written_by_print_XSF\n");
  fprintf(fxsf," DATAGRID_3D_DENSITY\n");
  fprintf(fxsf, "\t\t%d\t%d\t%d\n", NGX, NGY, NGZ);
  /*shift grid is set to 0 0 0 */
  fprintf(fxsf, "%lf\t%lf\t%lf\n", 0.0, 0.0, 0.0); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_ax, UC->ang_ay, UC->ang_az); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_bx, UC->ang_by, UC->ang_bz); 
  fprintf(fxsf, "\t%.10f\t%.10f\t%.10f\n", UC->ang_cx, UC->ang_cy, UC->ang_cz); 
  line_counter=0; 
  for(jz=0;jz<(NGZ);jz++) {
    for(jy=0;jy<(NGY);jy++) {
      for(jx=0;jx<(NGX);jx++) {
        line_counter++;
        fprintf(fxsf, "\t%.10lf" , BIN->real_grid[jx][jy][jz]);
        if(line_counter==6) {
          fprintf(fxsf, "\n");
          line_counter=0;
        }
      }
    }
  }
  fprintf(fxsf, "\nEND_DATAGRID_3D\n" );
  fprintf(fxsf, "END_BLOCK_DATAGRID3D" );
  fclose(fxsf);
  
  /*free allocated vars*/
  Xcart = FreeMemory_oneD_double(Xcart);
  Ycart = FreeMemory_oneD_double(Ycart);
  Zcart = FreeMemory_oneD_double(Zcart);
  BIN->real_grid = FreeMemory_threeD_double(BIN->real_grid, NGX, NGY);
} 
/*END of outputXSF function*/


void modify_line(char *line, Modification modifications[], int num_modifications)
{
  int max_line = 1024;
  int i;
  char modified_line[max_line];

  /*search through number of modifications to be made*/
  for (i=0;i<num_modifications;i++) {
    /*search for the word that should be modified*/
    if (strstr(line, modifications[i].search) != NULL) {
      /*if found, make the modification - copy the line with the modification before*/
      snprintf(modified_line, sizeof(modified_line), "%s%s", modifications[i].replace, line);
      /*replace the original line with the modified version*/
      strcpy(line, modified_line);
      break;
    }
  } 
}
  
void modify_abinitin(char in_filename[100], char out_filename[100])
{
  const int num_modifications = 8;
  FILE* abin;
  FILE* mabin;
  int max_line = 1024;
  char line[max_line];

  /*initialize modifications structure*/
  Modification modifications[] = {
{"toldfe", "!"},
{"occopt", "!"},
{"istwfk", "!"},
{"nshiftk", "!"},
{"shiftk", "!"},
{"ngkpt", "!"},
{"usekden", "!"},
{"prtkden", "!"},
  };

  /*open existing abinit in file to read*/
  printf("Modifying %s\n.", in_filename);
  abin=fopen(in_filename, "r");
  if (abin==NULL) {
    printf("%s not found. \n", in_filename);
    exit(0);
  }
  /*open a new file to modify the in file*/
  mabin=fopen(out_filename, "w");
  if (mabin==NULL) {
    printf("%s not found. \n", out_filename);
    exit(0);
  }

  while (fgets(line, sizeof(line), abin) != NULL) {
    modify_line(line, modifications, num_modifications);
    fputs(line, mabin);
  }

  fclose(abin); 
  fclose(mabin); 
      
} //END of  modify_abinitin function
    

