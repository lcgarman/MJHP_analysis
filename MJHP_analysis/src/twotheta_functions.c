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
#include "twotheta_functions.h"

void calculate_powder_pattern(FILE *flog, TwoTheta * TTH, BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC, Symmetry* SYM) 
{
  int ngfftx, ngffty, ngfftz;
  int h0, k0, l0;
  int h1, k1, l1;
  double F_hkl_val;
  double kx, ky, kz;
  double mag_G;
  double d_hkl_val;
  double lambda=0.71073; /*wavelenght of Mo radiation in ang*/
  double rad_theta;
  double deg_theta;
  double twotheta_val;
  int sym;
  int nsymor;
  int match;
  int j;
  int symm_nhkl;
  double F_diff;
  int ngfft_size;
  int * rflc_mult;
  int * H_arr;
  int * K_arr;
  int * L_arr;
  int H_symm, K_symm, L_symm;
  double tol;
  double * two_theta;
  double * d_hkl;
  double * F_hkl;
  int* hkl_sum;
  int hkl_sum_new;

  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  nsymor = SYM->nsymor;
  symm_nhkl = 0;
  tol = 1e-10;
  
  /*Initialize arrays*/
  rflc_mult = NULL; 
  H_arr = NULL; 
  K_arr = NULL;
  L_arr = NULL;
  hkl_sum = NULL;
  two_theta = NULL;
  d_hkl = NULL;
  F_hkl = NULL;
  /*allocate memory for temporary arrays*/
  ngfft_size = (ngfftx*ngffty*ngfftz);
  rflc_mult = AllocateMemory_oneD_int(rflc_mult, ngfft_size);
  H_arr = AllocateMemory_oneD_int(H_arr, ngfft_size);
  K_arr = AllocateMemory_oneD_int(K_arr, ngfft_size);
  L_arr = AllocateMemory_oneD_int(L_arr, ngfft_size);
  hkl_sum = AllocateMemory_oneD_int(hkl_sum, ngfft_size);
  two_theta = AllocateMemory_oneD_double(two_theta, ngfft_size);
  d_hkl = AllocateMemory_oneD_double(d_hkl, ngfft_size);
  F_hkl = AllocateMemory_oneD_double(F_hkl, ngfft_size);
  /*end of allocation*/

  /*Begin finding powder pattern*/
  fprintf(flog,  "\nBegin Calculating Powder Pattern\n");

  for (h0=0;h0<ngfftx;h0++) {
	for (k0=0;k0<ngffty;k0++) {
	  for (l0=0;l0<ngfftz;l0++) {
		if ((h0==0)&&(k0==0)&&(l0==0)) continue;

        F_hkl_val = gsl_complex_abs2(BIN->rec_grid[h0][k0][l0]);
        /*if the structure factor is zero skip*/
		if (F_hkl_val == 0.0) continue;

        /*assign actual HKL reflection values*/
		h1 = GRD->h_grid[h0][k0][l0];
		k1 = GRD->k_grid[h0][k0][l0];
		l1 = GRD->l_grid[h0][k0][l0];

        /*find the twotheta value for each reflection*/
		kx = h1*UC->ang_ax_star+k1*UC->ang_bx_star+l1*UC->ang_cx_star;
		ky = h1*UC->ang_ay_star+k1*UC->ang_by_star+l1*UC->ang_cy_star;
		kz = h1*UC->ang_az_star+k1*UC->ang_bz_star+l1*UC->ang_cz_star;
		mag_G = sqrt(kx*kx + ky*ky + kz*kz);
		d_hkl_val = (2.0*PI)/mag_G;
		if (d_hkl_val < 0.5*lambda) continue;
		rad_theta = asin((lambda)/(2*d_hkl_val));
		deg_theta = rad_theta * 180.0/PI;
		twotheta_val = deg_theta*2.0;
        /*if twotheta is greater than 60 skip*/
		if (twotheta_val > 60) continue;

        /*find symmetric hkl reflections*/
        match = 0;
		fprintf(flog,  "HKL = %d %d %d\n", h1, k1, l1);
        for(sym=0;sym<nsymor;sym++) {
          if (match == 1) break;
	      H_symm = SYM->symor[0][0][sym]*h1 + SYM->symor[1][0][sym]*k1 + SYM->symor[2][0][sym]*l1;
	      K_symm = SYM->symor[0][1][sym]*h1 + SYM->symor[1][1][sym]*k1 + SYM->symor[2][1][sym]*l1;
    	  L_symm = SYM->symor[0][2][sym]*h1 + SYM->symor[1][2][sym]*k1 + SYM->symor[2][2][sym]*l1;
	      for (j=0;j<symm_nhkl;j++) {
		    F_diff = F_hkl_val - F_hkl[j];
            /*search for symmetry related kx, ky, kz*/
            if ((H_symm==H_arr[j])&&(K_symm==K_arr[j])&&(L_symm==L_arr[j])&&(fabs(F_diff)<tol)) {
              rflc_mult[j]++;
              match = 1;
              fprintf(flog,  "\tMATCH(%e): %d %d %d F_hkl=%e\t == %d %d %d F_hkl=%e\n",F_diff, h1, k1, l1, F_hkl_val, H_arr[j], K_arr[j], L_arr[j], F_hkl[j]);
              hkl_sum_new = h1+k1+l1;
              if (hkl_sum_new>hkl_sum[j]) {
                H_arr[j] = h1;
                K_arr[j] = k1;
                L_arr[j] = l1;
                hkl_sum[j] = hkl_sum_new;
			  }
              break;
            }
          }
        }
        if (match == 0) {
	      H_arr[symm_nhkl] = h1;
	      K_arr[symm_nhkl] = k1;
		  L_arr[symm_nhkl] = l1;
		  two_theta[symm_nhkl] = twotheta_val;
		  d_hkl[symm_nhkl] = d_hkl_val;
		  F_hkl[symm_nhkl] = F_hkl_val;
		  rflc_mult[symm_nhkl] = 1;
          hkl_sum[symm_nhkl] = h1+k1+l1;
		  symm_nhkl++;
		}
	  }  //END l loop
	}  //END k loop
  }  //END h loop
  
  /*allocate memory for cpp arrays*/
  TTH->nrflc = symm_nhkl;
  TTH->rflc_mult = AllocateMemory_oneD_int(TTH->rflc_mult, TTH->nrflc);
  TTH->rflc_H = AllocateMemory_oneD_int(TTH->rflc_H, TTH->nrflc);
  TTH->rflc_K = AllocateMemory_oneD_int(TTH->rflc_K, TTH->nrflc);
  TTH->rflc_L = AllocateMemory_oneD_int(TTH->rflc_L, TTH->nrflc);
  TTH->two_theta = AllocateMemory_oneD_double(TTH->two_theta, TTH->nrflc);
  TTH->d_hkl = AllocateMemory_oneD_double(TTH->d_hkl, TTH->nrflc);
  TTH->F_hkl = AllocateMemory_oneD_double(TTH->F_hkl, TTH->nrflc);

  fprintf(flog,  "\nStoring Reflections:\n");
  /*store vaiables in PXRD struct*/
  for (j=0;j<symm_nhkl;j++) {
    TTH->rflc_H[j] = H_arr[j];
    TTH->rflc_K[j] = K_arr[j];
    TTH->rflc_L[j] = L_arr[j];
    TTH->two_theta[j] = two_theta[j];
    TTH->d_hkl[j] = d_hkl[j];
    TTH->F_hkl[j] = F_hkl[j];
    TTH->rflc_mult[j] = rflc_mult[j];
    fprintf(flog,  "\t%d HKL = %d %d %d\t mult=%d\t F_hkl = %lf\n", j, H_arr[j], K_arr[j], L_arr[j], rflc_mult[j], F_hkl[j]); 
  }

  /*free temporary allocated variables*/
  rflc_mult = FreeMemory_oneD_int(rflc_mult);
  H_arr = FreeMemory_oneD_int(H_arr);
  K_arr = FreeMemory_oneD_int(K_arr);
  L_arr = FreeMemory_oneD_int(L_arr);
  hkl_sum = FreeMemory_oneD_int(hkl_sum);
  two_theta = FreeMemory_oneD_double(two_theta);
  d_hkl = FreeMemory_oneD_double(d_hkl);
  F_hkl = FreeMemory_oneD_double(F_hkl);

  /*Free HKL grids wont need again*/
  GRD->h_grid = FreeMemory_threeD_int(GRD->h_grid, ngfftx, ngffty);
  GRD->k_grid = FreeMemory_threeD_int(GRD->k_grid, ngfftx, ngffty);
  GRD->l_grid = FreeMemory_threeD_int(GRD->l_grid, ngfftx, ngffty);
  
} //END of symmetry_reduce_pxrd

void fold_reflections_toBZ(FILE * flog, TwoTheta * TTH) 
{
  int nrflc;
  int n;

  nrflc = TTH->nrflc;
  /*allocate memory for temporary variables*/
  TTH->hpw = AllocateMemory_oneD_int(TTH->hpw, nrflc);
  TTH->kpw = AllocateMemory_oneD_int(TTH->kpw, nrflc);
  TTH->lpw = AllocateMemory_oneD_int(TTH->lpw, nrflc);
  TTH->BZkpt = AllocateMemory_twoD_double(TTH->BZkpt, nrflc, 3);
  /*end of allocation and initilization*/

  fprintf(flog,  "\nFolding PXRD Reflections Back into Brillouin Zone.\n");  

  for (n=0;n<nrflc;n++) {
    TTH->BZkpt[n][0] = (double) TTH->rflc_H[n]/2.0;
    TTH->BZkpt[n][1] = (double) TTH->rflc_K[n]/2.0;
    TTH->BZkpt[n][2] = (double) TTH->rflc_L[n]/2.0;
    TTH->hpw[n] = 0;
    TTH->kpw[n] = 0;
    TTH->lpw[n] = 0;
    while (TTH->BZkpt[n][0] > 0.5) {
      TTH->hpw[n]++;
      TTH->BZkpt[n][0] = TTH->BZkpt[n][0] - 1;
    }
    while (TTH->BZkpt[n][0] < -0.5) {
      TTH->hpw[n]--;
      TTH->BZkpt[n][0] = TTH->BZkpt[n][0] + 1;
    }
    while (TTH->BZkpt[n][1] > 0.5) {
      TTH->kpw[n]++;
      TTH->BZkpt[n][1] = TTH->BZkpt[n][1] - 1;
    }
    while (TTH->BZkpt[n][1] < -0.5) {
      TTH->kpw[n]--;
      TTH->BZkpt[n][1] = TTH->BZkpt[n][1] + 1;
    }
    while (TTH->BZkpt[n][2] > 0.5) {
      TTH->lpw[n]++;
      TTH->BZkpt[n][2] = TTH->BZkpt[n][2] - 1;
    }
    while (TTH->BZkpt[n][2] < -0.5) {
      TTH->lpw[n]--;
      TTH->BZkpt[n][2] = TTH->BZkpt[n][2] + 1;
    }
  }

} //END of fold_reflection_toBZ 

void symmetry_folded_reflections(FILE * flog, char filename[100], TwoTheta * TTH, Symmetry * SYM) 
{
  int nrflc;
  int nsymor;
  int sym;
  int match;
  int* kpt_match;
  double X_symm, Y_symm, Z_symm;
  int nsym_bzk;
  int n, i;
  double* kx_sym;
  double* ky_sym;
  double* kz_sym;
  double pw_x, pw_y, pw_z;
  int k;
  double kx, ky, kz;
  FILE* mabin;

  nrflc = TTH->nrflc;
  nsymor = SYM->nsymor;

  /*Initialize local arrays*/
  kx_sym = NULL;
  ky_sym = NULL;
  kz_sym = NULL;
  kpt_match = NULL;
  /*Dynamically Allocate memory for arrays in PXRD*/
  kx_sym = AllocateMemory_oneD_double(kx_sym, nrflc);
  ky_sym = AllocateMemory_oneD_double(ky_sym, nrflc);
  kz_sym = AllocateMemory_oneD_double(kz_sym, nrflc);
  kpt_match = AllocateMemory_oneD_int(kpt_match, nrflc);
  TTH->BZkpt_sym = AllocateMemory_twoD_double(TTH->BZkpt_sym, nrflc, 3);
  TTH->hpw_sym = AllocateMemory_oneD_int(TTH->hpw_sym, nrflc);
  TTH->kpw_sym = AllocateMemory_oneD_int(TTH->kpw_sym, nrflc);
  TTH->lpw_sym = AllocateMemory_oneD_int(TTH->lpw_sym, nrflc);
  TTH->rflc_H_sym = AllocateMemory_oneD_int(TTH->rflc_H_sym, nrflc);
  TTH->rflc_K_sym = AllocateMemory_oneD_int(TTH->rflc_K_sym, nrflc);
  TTH->rflc_L_sym = AllocateMemory_oneD_int(TTH->rflc_L_sym, nrflc);
  /*end of allocation*/

  fprintf(flog,  "\nApplying Symmetry to minimize kpts needed.\n");
  fprintf(flog, "\nApplying Symmetry to Folded Reflections.\n");
  /*intializing symmetric arrays*/
  kx_sym[0] = TTH->BZkpt[0][0];
  ky_sym[0] = TTH->BZkpt[0][1];
  kz_sym[0] = TTH->BZkpt[0][2];
  nsym_bzk = 1;

  /*Begin Finding symmetry related kpts*/
  for (n=0;n<nrflc;n++) {
    match=0;
    for (sym=0;sym<nsymor;sym++) {
      if (match==1) break;
	  X_symm = (double) SYM->symor[0][0][sym]*TTH->BZkpt[n][0] + (double) SYM->symor[1][0][sym]*TTH->BZkpt[n][1] + (double) SYM->symor[2][0][sym]*TTH->BZkpt[n][2];
	  Y_symm = (double) SYM->symor[0][1][sym]*TTH->BZkpt[n][0] + (double) SYM->symor[1][1][sym]*TTH->BZkpt[n][1] + (double) SYM->symor[2][1][sym]*TTH->BZkpt[n][2];
	  Z_symm = (double) SYM->symor[0][2][sym]*TTH->BZkpt[n][0] + (double) SYM->symor[1][2][sym]*TTH->BZkpt[n][1] + (double) SYM->symor[2][2][sym]*TTH->BZkpt[n][2];
	  //printf("\t%lf %lf %lf\t -> \t %lf %lf %lf\n", TTH->BZkpt[np1][0], TTH->BZkpt[np1][1],TTH->BZkpt[np1][2], X_symm, Y_symm, Z_symm);
	  
      for (i=0;i<nsym_bzk;i++) {
		if ((X_symm==kx_sym[i])&&(Y_symm==ky_sym[i])&&(Z_symm==kz_sym[i])) {
          match = 1;
		  TTH->rflc_H_sym[n] = SYM->symor[0][0][sym]*TTH->rflc_H[n] + SYM->symor[1][0][sym]*TTH->rflc_K[n] + SYM->symor[2][0][sym]*TTH->rflc_L[n];
		  TTH->rflc_K_sym[n] = SYM->symor[0][1][sym]*TTH->rflc_H[n] + SYM->symor[1][1][sym]*TTH->rflc_K[n] + SYM->symor[2][1][sym]*TTH->rflc_L[n];
		  TTH->rflc_L_sym[n] = SYM->symor[0][2][sym]*TTH->rflc_H[n] + SYM->symor[1][2][sym]*TTH->rflc_K[n] + SYM->symor[2][2][sym]*TTH->rflc_L[n];
          kpt_match[n] = i;
          break;
        }
      }
      if (((sym)==(nsymor-1))&&(match==0)) {
		kx_sym[nsym_bzk] = TTH->BZkpt[n][0];
		ky_sym[nsym_bzk] = TTH->BZkpt[n][1];
		kz_sym[nsym_bzk] = TTH->BZkpt[n][2];
		TTH->rflc_H_sym[n] = TTH->rflc_H[n];
		TTH->rflc_K_sym[n] = TTH->rflc_K[n];
		TTH->rflc_L_sym[n] = TTH->rflc_L[n];
        kpt_match[n] = nsym_bzk;
		nsym_bzk++;
      }
	}
  }

  fprintf(flog,  "\t Symmetrized kpts and Corresponding Reflections:\n");
  for (n=0;n<nrflc;n++) {
    pw_x = (double) TTH->rflc_H_sym[n]/2.0;
    pw_y = (double) TTH->rflc_K_sym[n]/2.0;
    pw_z = (double) TTH->rflc_L_sym[n]/2.0;
    k = kpt_match[n];
    kx = kx_sym[k];
    ky = ky_sym[k];
    kz = kz_sym[k];
    TTH->BZkpt_sym[n][0] = kx;
    TTH->BZkpt_sym[n][1] = ky;
    TTH->BZkpt_sym[n][2] = kz;
    TTH->hpw_sym[n] = pw_x - kx; 
    TTH->kpw_sym[n] = pw_y - ky; 
    TTH->lpw_sym[n] = pw_z - kz; 
	fprintf(flog,  "\t %lf %lf %lf\t%lf %lf %lf %d %d %d\n", pw_x, pw_y, pw_z, TTH->BZkpt_sym[n][0], TTH->BZkpt_sym[n][1], TTH->BZkpt_sym[n][2], TTH->hpw_sym[n], TTH->kpw_sym[n], TTH->lpw_sym[n]);
  }

  /*print into modified abinit in file*/
  mabin=fopen(filename, "a");
  if (mabin==NULL) {
    fprintf(flog, "%s not found. \n", filename);
    exit(0);
  }
  fprintf(mabin,"\n!---- MJHP 2THETA VARS ----!\n");
  fprintf(mabin,"      istwfk %d*1\n", nsym_bzk);
  fprintf(mabin,"      tolwfr 0.0000000000001\n");
  fprintf(mabin,"      iscf -2 \n");
  fprintf(mabin,"      kptopt 0\n");
  fprintf(mabin,"      nkpt %d\n", nsym_bzk);
  fprintf(mabin,"      kpt\n");
  for (i=0;i<nsym_bzk;i++) {
    fprintf(mabin,"      %lf %lf %lf\n", kx_sym[i], ky_sym[i], kz_sym[i]);
  }  
  fclose(mabin);

  /*Free local memory*/
  kx_sym = FreeMemory_oneD_double(kx_sym);
  ky_sym = FreeMemory_oneD_double(ky_sym);
  kz_sym = FreeMemory_oneD_double(kz_sym);
  kpt_match = FreeMemory_oneD_int(kpt_match);

} //END of Symmetry_pxrd

void print_reflections(char filename[200], TwoTheta * TTH)
{
  FILE * frflc;
  int n;
  int nrflc;
  double intensity;
  nrflc = TTH->nrflc;
  
  /*Printf Symmetrized kpts*/
  printf("\nPrinting Reflection Information into %s\n", filename);
  frflc=fopen(filename,"w"); /*open the filename in append mode*/
  if(frflc==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  fprintf(frflc, "nrflc %d\n", nrflc);
  /*print powder pattern */
  for (n=0;n<nrflc;n++) {
    intensity = TTH->F_hkl[n] * TTH->rflc_mult[n];
    fprintf(frflc, "%d\t %lf %lf %d %d %d %d %lf %lf\n", n, TTH->d_hkl[n], TTH->two_theta[n], TTH->rflc_H[n], TTH->rflc_K[n], TTH->rflc_L[n], TTH->rflc_mult[n], TTH->F_hkl[n], intensity);
  }
  
  /*printf symmetrized and folded BZkpts and HKLpw*/
  fprintf(frflc, "\n");
  fprintf(frflc, "#__H_K_L___kx_ky_kx____hpw_kpw_lpw\n");
  for  (n=0;n<nrflc;n++) {
    fprintf(frflc, "%d\t%d %d %d\t%lf %lf %lf\t%d %d %d\n", n, TTH->rflc_H_sym[n], TTH->rflc_K_sym[n], TTH->rflc_L_sym[n], TTH->BZkpt_sym[n][0], TTH->BZkpt_sym[n][1], TTH->BZkpt_sym[n][2], TTH->hpw_sym[n], TTH->kpw_sym[n], TTH->lpw_sym[n]);
  }

  fclose(frflc);
  
  /*Free allocated variables*/
  TTH->rflc_mult = FreeMemory_oneD_int(TTH->rflc_mult);
  TTH->rflc_H = FreeMemory_oneD_int(TTH->rflc_H);
  TTH->rflc_K = FreeMemory_oneD_int(TTH->rflc_K);
  TTH->rflc_L = FreeMemory_oneD_int(TTH->rflc_L);
  TTH->two_theta = FreeMemory_oneD_double(TTH->two_theta);
  TTH->d_hkl = FreeMemory_oneD_double(TTH->d_hkl);
  TTH->F_hkl = FreeMemory_oneD_double(TTH->F_hkl);
  TTH->hpw = FreeMemory_oneD_int(TTH->hpw);
  TTH->kpw = FreeMemory_oneD_int(TTH->kpw);
  TTH->lpw = FreeMemory_oneD_int(TTH->lpw);
  TTH->BZkpt = FreeMemory_twoD_double(TTH->BZkpt, nrflc);
  TTH->BZkpt_sym = FreeMemory_twoD_double(TTH->BZkpt_sym, nrflc);
  TTH->hpw_sym = FreeMemory_oneD_int(TTH->hpw_sym);
  TTH->kpw_sym = FreeMemory_oneD_int(TTH->kpw_sym);
  TTH->lpw_sym = FreeMemory_oneD_int(TTH->lpw_sym);
  TTH->rflc_H_sym = FreeMemory_oneD_int(TTH->rflc_H_sym);
  TTH->rflc_K_sym = FreeMemory_oneD_int(TTH->rflc_K_sym);
  TTH->rflc_L_sym = FreeMemory_oneD_int(TTH->rflc_L_sym);
} //END of Append_PXRDfile

void read_reflections(char filename[200], TwoTheta * TTH)
{
  FILE * frflc;
  int n;
  int nrflc;
  int rflc_no;
  double d_hklv;
  double two_thetav;
  int rflc_Hv, rflc_Kv, rflc_Lv;
  int rflc_multv;
  double F_hklv;
  double intensity;
  char strv[50];
  double kx, ky, kz;
  int hpwv, kpwv, lpwv;
  
  /*Printf Symmetrized kpts*/
  printf("\nReading Reflection Information From: %s\n", filename);
  frflc=fopen(filename,"r"); /*open the filename in append mode*/
  if(frflc==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  fscanf(frflc, "%s", strv);
  fscanf(frflc, "%d", &nrflc);
  TTH->nrflc = nrflc;

  /*allocate memroy for reflection information*/
  TTH->rflc_mult = AllocateMemory_oneD_int(TTH->rflc_mult, nrflc);
  TTH->rflc_H = AllocateMemory_oneD_int(TTH->rflc_H, nrflc);
  TTH->rflc_K = AllocateMemory_oneD_int(TTH->rflc_K, nrflc);
  TTH->rflc_L = AllocateMemory_oneD_int(TTH->rflc_L, nrflc);
  TTH->two_theta = AllocateMemory_oneD_double(TTH->two_theta, nrflc);
  TTH->BZkpt = AllocateMemory_twoD_double(TTH->BZkpt, nrflc, 3);
  TTH->hpw = AllocateMemory_oneD_int(TTH->hpw, nrflc);
  TTH->kpw = AllocateMemory_oneD_int(TTH->kpw, nrflc);
  TTH->lpw = AllocateMemory_oneD_int(TTH->lpw, nrflc);
  /*end of allocation*/

  /*read powder pattern */
  for (n=0;n<nrflc;n++) {
    fscanf(frflc, "%d", &rflc_no);
    fscanf(frflc, "%lf", &d_hklv);
    fscanf(frflc, "%lf", &two_thetav);
    fscanf(frflc, "%d %d %d ", &rflc_Hv, &rflc_Kv, &rflc_Lv);
    fscanf(frflc, "%d", &rflc_multv);
    fscanf(frflc, "%lf", &F_hklv);
    fscanf(frflc, "%lf", &intensity);
    /*store necessary variables in TTH struct*/
    TTH->two_theta[n] = two_thetav;
    TTH->rflc_mult[n] = rflc_multv;
  }
  
  /*read symmetrized and folded BZkpts and HKLpw*/
  fscanf(frflc, "%s", strv);
  for  (n=0;n<nrflc;n++) {
    fscanf(frflc, "%d", &rflc_no);
    fscanf(frflc, "%d %d %d", &rflc_Hv, &rflc_Kv, &rflc_Lv);
    fscanf(frflc, "%lf %lf %lf", &kx, &ky, &kz);
    fscanf(frflc, "%d %d %d", &hpwv, &kpwv, &lpwv);
    TTH->rflc_H[n] = rflc_Hv;
    TTH->rflc_K[n] = rflc_Kv;
    TTH->rflc_L[n] = rflc_Lv;
    TTH->BZkpt[n][0] = kx;
    TTH->BZkpt[n][1] = ky;
    TTH->BZkpt[n][2] = kz;
    TTH->hpw[n] = hpwv;
    TTH->kpw[n] = kpwv;
    TTH->lpw[n] = lpwv;
  }

  fclose(frflc);
}//END of read in rflc file
  
void concatinate_twotheta_potential(EnergyContribution * ECON, EnergyStep * ESTP, TwoTheta *TTH) 
{
  int nEstep;
  int dE;
  int nrflc;
  int n;
  double reflection_total;

  nEstep = ESTP->nEstep;
  nrflc = TTH->nrflc;
  /*Allocate Memory*/
  ECON->rflc_total = AllocateMemory_twoD_double(ECON->rflc_total, nEstep, nrflc);
  for (dE=0;dE<nEstep;dE++) {
    for (n=0;n<nrflc;n++) {
      ECON->rflc_total[dE][n] = 0.0;
    }
  }
  
  printf("\nConcatinating local and nonlocal potential grids.\n");
  
  /*add local and nonlocal grids to find total potential energy*/
  for (n=0;n<nrflc;n++) {
    reflection_total = 0.0;
    for (dE=0;dE<nEstep;dE++) {
	  ECON->rflc_total[dE][n] = ECON->rflc_local[dE][n] + ECON->rflc_nonlocal[dE][n];
	  reflection_total += ECON->rflc_total[dE][n];
    }
    printf( "Total Potential Energy for rflc %d = %lf \n", n, reflection_total);
  }
  
  /*free memory for local and nonlocal*/
  ECON->rflc_local = FreeMemory_twoD_double(ECON->rflc_local, nEstep);
  ECON->rflc_nonlocal = FreeMemory_twoD_double(ECON->rflc_nonlocal, nEstep);
}

