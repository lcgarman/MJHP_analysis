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
#include "mjhpHKL_local_potential.h"

void mjhpHKL_local_potential(VectorIndices *VECT, NumberGrid *GRD, EnergyStep *ESTP, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, EnergyContribution *ECON)
{ 
  int nkpt;
  int kptno;
  int nband;
  int band;
  int npw;
  int nEstep;
  double bandE_min;
  double fermi;
  int Emesh;
  int ngfftx, ngffty, ngfftz;
  double bandE;
  int dE;
  double occ;
  int pw1, pw2;
  int h1, k1, l1;
  int h2, k2, l2;
  int delta_h, delta_k, delta_l;
  double pw1_x, pw1_y, pw1_z;
  double pw2_x, pw2_y, pw2_z;
  double ang_pw1x, ang_pw1y, ang_pw1z;
  double ang_pw2x, ang_pw2y, ang_pw2z;
  double mag_pw1, mag_pw2;
  double mag_diff;
  double exponent;
  double broad;
  double minr, maxr;
  int j;
  int nHKL;
  int hkl_match;
  int H_match, K_match, L_match;
  double wtk;
  double sigma;
  double one_wavecoef_RE, one_wavecoef_IM;
  double two_wavecoef_RE, two_wavecoef_IM;
  gsl_complex c1, c1_star;
  gsl_complex c2;
  gsl_complex c1star_c2;
  gsl_complex Vhkl;
  gsl_complex complexnums;
  double noncomplex;
  gsl_complex potentialE;
  double total_local;
  
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  nEstep = ESTP->nEstep;
  bandE_min = ESTP->bandE_min;
  fermi  = WFK->fermi;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  minr = VECT->minr;
  maxr = VECT->maxr;
  nHKL = VECT->nHKL;
  sigma = SIGMA;
  Emesh = EMESH;
  total_local = 0.0;
 
  /*Allocate memory for variables in this function*/
  ECON->local = AllocateMemory_oneD_double(ECON->local, nEstep);

  printf( "\nCalculating Local Potential Energy:\n");
  /*NOW start Calculating Potential Energy contribution for each kpt*/
  for (kptno=0;kptno<nkpt;kptno++){ 
    printf( "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

    /*loop over first planewaves in pair*/
	for(pw1=0;pw1<npw;pw1++) {
      /*reduced coordinates of planewave*/
	  h1 = WFK->kg[kptno][pw1][0];
	  k1 = WFK->kg[kptno][pw1][1];
	  l1 = WFK->kg[kptno][pw1][2];
	  /*Find magnitude of pw1 in inverse angstroms*/ 
	  pw1_x = WFK->kpt[0][kptno] + (double) h1;
	  pw1_y = WFK->kpt[1][kptno] + (double) k1;
	  pw1_z = WFK->kpt[2][kptno] + (double) l1;
	  ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
	  ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
	  ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
	  mag_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
      /*if the magnitude of pw1 is outside of shell go to next iteration*/
      if ((mag_pw1>maxr)||(mag_pw1<minr)) continue;
	  
      /*now loop over second planewave in pair*/
	  for(pw2=0;pw2<npw;pw2++) {
		h2 = WFK->kg[kptno][pw2][0];
		k2 = WFK->kg[kptno][pw2][1];
		l2 = WFK->kg[kptno][pw2][2];
	    /*Find magnitude of pw2*/ 
		pw2_x = WFK->kpt[0][kptno] + (double) h2;
		pw2_y = WFK->kpt[1][kptno] + (double) k2;
		pw2_z = WFK->kpt[2][kptno] + (double) l2;
		ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
	    mag_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
        /*if the magnitude of pw2 is outside of shell go to next iteration*/
        if ((mag_pw2>maxr)||(mag_pw2<minr)) continue;
		
        /*if pw1 and pw2 both lie in shell; find the difference bw them*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;
        hkl_match = 0;
        /*check all symmetry equialent HKL to see if pw1-pw2 matches these indices*/
        for(j=0;j<nHKL;j++) {
		  if ((delta_h==VECT->H_arr[j])&&(delta_k==VECT->K_arr[j])&&(delta_l==VECT->L_arr[j])) {
            hkl_match = 1;
            H_match = VECT->H_arr[j];
            K_match = VECT->K_arr[j];
            L_match = VECT->L_arr[j];
            break; /*if we find a match stop looking; exit loop*/
		  }
        }
        /*if no match is found go to next iteration*/
        if (hkl_match == 0) continue;
         
        /*print out pw matches and coordinates*/
        printf( "\t H K L = %d %d %d\n", H_match, K_match, L_match);
        printf( "\t\t pw1 = %lf %lf %lf |pw1|=%lf\n", pw1_x, pw1_y, pw1_z, mag_pw1);
        printf( "\t\t pw2 = %lf %lf %lf |pw2|=%lf\n", pw2_x, pw2_y, pw2_z, mag_pw2);
 
        /*make delta_h positive*/
		if (delta_h < 0) delta_h+=ngfftx;
		if (delta_k < 0) delta_k+=ngffty;
		if (delta_l < 0) delta_l+=ngfftz;

        /*calculate broadening*/
        mag_diff = mag_pw1 - mag_pw2;
        exponent = -(mag_diff*mag_diff)/sigma;
        broad = exp(exponent);
        printf( "\t\tbroadening = %lf\n", broad);

		/*loop over bands to find wavefunction coeff and dE*/
		for(band=0;band<nband;band++) {
		  /*Determine Energy range to store contribution*/
		  bandE = (((WFK->eigen[kptno][band]-fermi)*HATOEV-bandE_min)*Emesh)+0.5;
		  dE = floor(bandE);
		  if (dE >= nEstep) continue; /*go to next iteration if outside of Erange*/ 
		  
		  /*Setting Band Occupation and Kpt Weight*/
		  //occ = WFK->occ[kptno][band];
		  occ = 1;
          wtk = WFK->wtk[kptno];
		  
		  /*Calculate potential energy for HKL from pw1 and pw2*/
		  one_wavecoef_RE = WFK->cg[kptno][band][pw1][1];
		  one_wavecoef_IM = WFK->cg[kptno][band][pw1][0];
		  c1 = gsl_complex_rect(one_wavecoef_RE, one_wavecoef_IM);
		  c1_star = gsl_complex_conjugate(c1);
		  two_wavecoef_RE = WFK->cg[kptno][band][pw2][1];
		  two_wavecoef_IM = WFK->cg[kptno][band][pw2][0];
		  c2 = gsl_complex_rect(two_wavecoef_RE, two_wavecoef_IM);
		  
		  c1star_c2 = gsl_complex_mul(c1_star, c2);
		  Vhkl = BIN->rec_grid[delta_h][delta_k][delta_l];
		  complexnums = gsl_complex_mul(c1star_c2, Vhkl);
		  noncomplex = occ*wtk*broad/UC->bohr_cellV;
		  potentialE = gsl_complex_mul_real(complexnums, noncomplex);

          /*store potential energy contributions*/
          ECON->local[dE] += GSL_REAL(potentialE);
		  total_local += GSL_REAL(potentialE);
		} /*END band->nband loop*/
	  } /*END pw2->npw loop*/
	} /*END pw1->npw loop*/
	printf( "   kpt %d\t local potential energy = %lf\n", kptno, total_local);
  } /*END: kpt loop*/
  printf("\tTotal Potential Energy = %lf\n", total_local);
  printf( "Total Potential Energy = %lf\n", total_local);


}   //END of mjhp_hkl_localpot function
