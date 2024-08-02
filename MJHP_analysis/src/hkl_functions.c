#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "hkl_functions.h"

void transform_HKL(MottJonesConditions * MJC) 
{
  char lattice[10];
  int len_lattice;
  char center;
  char primitive = 'P';
  char bodycenter = 'I';
  char facecenter = 'F';
  char acenter = 'A';
  char bcenter = 'B';
  char ccenter = 'C';
  int H_symm, K_symm, L_symm;
  int H, K, L;
  
  H = MJC->nH;
  K = MJC->nK;
  L = MJC->nL;
  /*separating lattice variable. ie - "cP" -> "c" and "P"*/
  strcpy(lattice, MJC->lattice);
  len_lattice = strlen(lattice);
  /*check lattice is right length*/
  if (len_lattice > 2) printf("ERROR: length of lattice (%d) is too long.\n", len_lattice);
  center = lattice[1];

  /*Transform to primitive cell*/
  if (center==primitive) {
    printf("Already in Primitive Centering\n");
  }
  else if (center==bodycenter) { 
    printf("Transforming HKL: I->P\n");
	H_symm = (-0.5*H)+(0.5*K)+(0.5*L);
	K_symm = (0.5*H)+(-0.5*K)+(0.5*L);
	L_symm = (0.5*H)+(0.5*K)+(-0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==facecenter) { 
    printf("Transforming HKL: F->P\n");
	H_symm = (0.0*H)+(0.5*K)+(0.5*L);
	K_symm = (0.5*H)+(0.0*K)+(0.5*L);
	L_symm = (0.5*H)+(0.5*K)+(0.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==acenter) { 
    printf("Transforming HKL: A->P\n");
	H_symm = (1.0*H)+(0.0*K)+(0.0*L);
	K_symm = (0.0*H)+(0.5*K)+(-0.5*L);
	L_symm = (0.0*H)+(0.5*K)+(0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==bcenter) { 
    printf("Transforming HKL: B->P\n");
	H_symm = (0.5*H)+(0.0*K)+(-0.5*L);
	K_symm = (0.0*H)+(1.0*K)+(0.0*L);
	L_symm = (0.5*H)+(0.0*K)+(0.5*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==ccenter) { 
    printf("Transforming HKL: C->P\n");
	H_symm = (0.5*H)+(-0.5*K)+(0.0*L);
	K_symm = (0.5*H)+(0.5*K)+(0.0*L);
	L_symm = (0.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else {
    printf("ERROR: No Code written for %s\n", MJC->lattice);
    printf("\t case sensitive (xY) \n");
    exit(0);
  }

  MJC->nH = H;
  MJC->nK = K;
  MJC->nL = L;
} //END of Transform_hklreflection_

void find_symmetric_hkl(VectorIndices* VECT, Symmetry* SYM, NumberGrid* GRD) 
{
  int j;
  int sym;
  int nsym;
  int H, K, L;
  int nHKL;
  int nomatch;
  int H_symm, K_symm, L_symm;
  int* H_arr;
  int* K_arr;
  int* L_arr;
  int ngfftx, ngffty, ngfftz;
  
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;
  nsym = SYM->nsym;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  H_arr = NULL;
  K_arr = NULL;
  L_arr = NULL;
  
  /*Allocate Memory for HKL arrays*/
  H_arr = AllocateMemory_oneD_int(H_arr, nsym);
  K_arr = AllocateMemory_oneD_int(K_arr, nsym);
  L_arr = AllocateMemory_oneD_int(L_arr, nsym);
  /*End of allocate memory*/

  printf( "\nSymmetrically equivalent HKL\n");
  /*assign first element of HKL array to input HKL*/
  H_arr[0] = H;
  K_arr[0] = K;
  L_arr[0] = L;
  nHKL = 1;
  /*Find symmetry equivalent reflections based off of nonsymmorphic symm elements*/
  for(sym=0;sym<nsym;sym++) {
	H_symm = SYM->symrel[0][0][sym]*H + SYM->symrel[1][0][sym]*K + SYM->symrel[2][0][sym]*L;
	K_symm = SYM->symrel[0][1][sym]*H + SYM->symrel[1][1][sym]*K + SYM->symrel[2][1][sym]*L;
	L_symm = SYM->symrel[0][2][sym]*H + SYM->symrel[1][2][sym]*K + SYM->symrel[2][2][sym]*L;
	nomatch = 0;
	//printf("no%d\t H = %d\t K = %d\t L = %d\n", sym, H_symm, K_symm, L_symm);
	//printf("\t %d*%d + %d*%d + %d*%d = %d\n",  SYM->symrel[0][0][sym], H,  SYM->symrel[1][0][sym], K, SYM->symrel[2][0][sym], L, H_symm); 
	//printf("\t %d*%d + %d*%d + %d*%d = %d\n",  SYM->symrel[0][1][sym], H,  SYM->symrel[1][1][sym], K, SYM->symrel[2][1][sym], L, K_symm); 
	//printf("\t %d*%d + %d*%d + %d*%d = %d\n",  SYM->symrel[0][2][sym], H,  SYM->symrel[1][2][sym], K, SYM->symrel[2][2][sym], L, L_symm); 

	/*search elements in array if you find a match go to next symm element; if not store*/
    for (j=0;j<nHKL;j++) {
	  if ((H_symm == H_arr[j])&&(K_symm == K_arr[j])&&(L_symm == L_arr[j])) {
	    break;
	  }
	  else nomatch++;
    }
	if (nomatch == nHKL) {
	  H_arr[nHKL] = H_symm;
	  K_arr[nHKL] = K_symm;
	  L_arr[nHKL] = L_symm;
	  nHKL++;
	}
  }
  
  printf( "\tMultiplicity of HKL Reflection = %d\n", nHKL);
  /*allocate memory for RFLC arrays*/
  VECT->H_arr = AllocateMemory_oneD_int(VECT->H_arr, nHKL);
  VECT->K_arr = AllocateMemory_oneD_int(VECT->K_arr, nHKL);
  VECT->L_arr = AllocateMemory_oneD_int(VECT->L_arr, nHKL);
  VECT->H_posarr = AllocateMemory_oneD_int(VECT->H_posarr, nHKL);
  VECT->K_posarr = AllocateMemory_oneD_int(VECT->K_posarr, nHKL);
  VECT->L_posarr = AllocateMemory_oneD_int(VECT->L_posarr, nHKL);
  /*end of allocation*/

  /*store variables in VectorIndices structure*/
  for (j=0;j<nHKL;j++) {
    VECT->H_arr[j] = H_arr[j];	
    if (H_arr[j] < 0) {
      VECT->H_posarr[j] = H_arr[j] + ngfftx; 
    }
    else {
      VECT->H_posarr[j] = H_arr[j]; 
    }
    VECT->K_arr[j] = K_arr[j]; 
    if (K_arr[j] < 0) {
      VECT->K_posarr[j] = K_arr[j] + ngffty; 
    }
    else {
      VECT->K_posarr[j] = K_arr[j]; 
    }
    VECT->L_arr[j] = L_arr[j];
    if (L_arr[j] < 0) {
      VECT->L_posarr[j] = L_arr[j] + ngfftz; 
    }
    else {
      VECT->L_posarr[j] = L_arr[j]; 
    }
	printf( "\t#%d:   %d(%d) %d(%d) %d(%d)\n", j+1, VECT->H_arr[j], VECT->H_posarr[j], VECT->K_arr[j], VECT->K_posarr[j], VECT->L_arr[j], VECT->L_posarr[j]);
  }
  VECT->nHKL = nHKL;

  /*free allocated variables*/
  H_arr = FreeMemory_oneD_int(H_arr);
  K_arr = FreeMemory_oneD_int(K_arr);
  L_arr = FreeMemory_oneD_int(L_arr);

} //END of Find_HKLsymmetry
  
void find_MJregion(VectorIndices *VECT, UnitCell *UC) 
{
  int H, K, L;
  double x0, y0, z0;
  double kx, ky, kz;
  double rad;
  double mag_Ghkl;
  double minr, maxr;
  
  printf( "\n Finding shell in reciprocal space to analyze.\n");
  /*Initialize Variables*/
  H = VECT->H;
  K = VECT->K;
  L = VECT->L;
  rad = 0.075;

  /*Find JZ face-center*/
  x0 = (double) H/2.0;
  y0 = (double) K/2.0;
  z0 = (double) L/2.0;
  /*calculate ang-1 coords of jzfc*/
  kx = x0*UC->ang_ax_star + y0*UC->ang_bx_star + z0*UC->ang_cx_star;
  ky = x0*UC->ang_ay_star + y0*UC->ang_by_star + z0*UC->ang_cy_star;
  kz = x0*UC->ang_az_star + y0*UC->ang_bz_star + z0*UC->ang_cz_star;
  /*calc magnitude/ distance of kx,ky,kx*/
  mag_Ghkl = sqrt(kx*kx+ky*ky+kz*kz);
  printf( "\tCenter of JZ face = %lf %lf %lf\t (%lf %lf %lf)ang-1 |G_hkl|=%lf ang-1\n", x0, y0, z0, kx, ky, kz, mag_Ghkl);
  
  minr = mag_Ghkl - rad;
  maxr = mag_Ghkl + rad;
  printf( "\tShell: inner radius = %lf ang-1\t outer radius = %lfang-1\n", minr, maxr);
  VECT->minr = minr;
  VECT->maxr = maxr;

} //END of find_MJregion 

void concatinate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP) 
{
  int nEstep;
  int dE;
  double total_potential;

  nEstep = ESTP->nEstep;
  /*Allocate Memory*/
  ECON->total_potential = AllocateMemory_oneD_double(ECON->total_potential, nEstep);
  for (dE=0;dE<nEstep;dE++) ECON->total_potential[dE] = 0.0;
  
  printf( "\nConcatinating local and nonlocal potential grids.\n");
  
  /*add local and nonlocal grids to find total potential energy*/
  total_potential = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    ECON->total_potential[dE] = ECON->local[dE] + ECON->nonlocal[dE];
    total_potential += ECON->total_potential[dE];
  }
  
  printf("\tTotal Potential Energy is %lf \n", total_potential);

  /*free memory for local and nonlocal*/
  ECON->local = FreeMemory_oneD_double(ECON->local);
  ECON->nonlocal = FreeMemory_oneD_double(ECON->nonlocal);
}

void integrate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP, AtomicVariables * ATM)
{
  int nEstep;
  int dE;
  double Eint_Ha;
  double Eint_eV_peratom;
  int dE_zero;
  int natom;
 
  nEstep = ESTP->nEstep;
  dE_zero = ESTP->dE_zero;
  natom = ATM->natom;

  printf("\nIntegrated the Total Potential Energy.\n");
  
  Eint_Ha = 0.0;
  for (dE=0;dE<nEstep;dE++) {
    if (dE > dE_zero) continue;
    Eint_Ha += ECON->total_potential[dE];
  }
  Eint_eV_peratom = (Eint_Ha*HATOEV)/((double)natom);

  printf( "\t iMJHP = %lf eV/atom\n", Eint_eV_peratom);

}

Conventional HKL_convert_toP(MottJonesConditions * MJC, int H, int K, int L)
{
  char lattice[10];
  int len_lattice;
  char center;
  char primitive = 'P';
  char bodycenter = 'I';
  char facecenter = 'F';
  char acenter = 'A';
  char bcenter = 'B';
  char ccenter = 'C';
  int H_symm, K_symm, L_symm;
  
  Conventional con;
  /*separating lattice variable. ie - "cP" -> "c" and "P"*/
  strcpy(lattice, MJC->lattice);
  len_lattice = strlen(lattice);
  /*check lattice is right length*/
  if (len_lattice > 2) printf("ERROR: length of lattice (%d) is too long.\n", len_lattice);
  center = lattice[1];

  /*Transform to primitive cell*/
  if (center==primitive) {
    printf("Already in Conventional Centering (P)\n");
  }
  else if (center==bodycenter) { 
    printf("Transforming HKL: P->I\n");
	H_symm = (0.0*H)+(1.0*K)+(1.0*L);
	K_symm = (1.0*H)+(0.0*K)+(1.0*L);
	L_symm = (1.0*H)+(1.0*K)+(0.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==facecenter) { 
    printf("Transforming HKL: P->F\n");
	H_symm = (-1.0*H)+(1.0*K)+(1.0*L);
	K_symm = (1.0*H)+(-1.0*K)+(1.0*L);
	L_symm = (1.0*H)+(1.0*K)+(-1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==acenter) { 
    printf("Transforming HKL: A->P\n");
	H_symm = (1.0*H)+(0.0*K)+(0.0*L);
	K_symm = (0.0*H)+(1.0*K)+(1.0*L);
	L_symm = (0.0*H)+(-1.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==bcenter) { 
    printf("Transforming HKL: P->B\n");
	H_symm = (1.0*H)+(0.0*K)+(1.0*L);
	K_symm = (0.0*H)+(1.0*K)+(0.0*L);
	L_symm = (-1.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else if (center==ccenter) { 
    printf("Transforming HKL: P->C\n");
	H_symm = (1.0*H)+(1.0*K)+(0.0*L);
	K_symm = (-1.0*H)+(1.0*K)+(0.0*L);
	L_symm = (0.0*H)+(0.0*K)+(1.0*L);
    H = H_symm;
    K = K_symm;
    L = L_symm;
  }
  else {
    printf("ERROR: No Code written for %s\n", MJC->lattice);
    printf("\t case sensitive (xY) \n");
    exit(0);
  }

  con.H = H;
  con.K = K;
  con.L = L;
  
  return(con);
} //END of primitive_to_conventional 

