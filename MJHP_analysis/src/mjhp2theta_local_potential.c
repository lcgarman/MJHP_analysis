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
#include "mjhp2theta_local_potential.h"

void mjhp2theta_local_potential(TwoTheta *TTH, EnergyStep *ESTP, NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, EnergyContribution *ECON)
{ /*Calculate the kinetic energy contribution by HKL value*/
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
  double mag_pw1, mag_pw2;
  int H, K, L;
  int hpw, kpw, lpw;
  int rflc_mult;
  double kx, ky, kz;
  int nrflc;
  int n;
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
  Emesh = EMESH;
  nrflc = TTH->nrflc;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
 
  /*Allocate memory for large variables in this function*/
  ECON->rflc_local = AllocateMemory_twoD_double(ECON->rflc_local, nEstep, nrflc);
  //for (dE=0;dE<nEstep;dE++) {
    //for (n=0;n<nrflc;n++) {
      //ECON->rflc_local[dE][n] = 0.0;
   // }
 // }
  total_local = 0.0;
  
  printf( "\nCalculating Local Potential Energy:\n");
  /*Enter loop of  reflections*/
  for (n=0;n<nrflc;n++) {
    kx = TTH->BZkpt[n][0];  
    ky = TTH->BZkpt[n][1];  
    kz = TTH->BZkpt[n][2];  
    hpw = TTH->hpw[n];
    kpw = TTH->kpw[n];
    lpw = TTH->lpw[n];
    H = TTH->rflc_H[n];
    K = TTH->rflc_K[n];
    L = TTH->rflc_L[n];
    rflc_mult = TTH->rflc_mult[n];
	/*NOW start Calculating local potential  Energy contribution*/
	for (kptno=0;kptno<nkpt;kptno++){ 
	  if ((kx!=WFK->kpt[0][kptno])||(ky!=WFK->kpt[1][kptno])||(kz!=WFK->kpt[2][kptno])) continue;
	  printf( "\treflection# %d\n", n);
	  npw = WFK->npw[kptno];
      /*loop over first planewave in pair*/
	  for(pw1=0;pw1<npw;pw1++) {
		h1 = WFK->kg[kptno][pw1][0];
		k1 = WFK->kg[kptno][pw1][1];
		l1 = WFK->kg[kptno][pw1][2];
		if ((h1!=hpw)||(k1!=kpw)||(l1!=lpw)) continue;
		pw1_x = WFK->kpt[0][kptno] + (double) h1;
		pw1_y = WFK->kpt[1][kptno] + (double) k1;
		pw1_z = WFK->kpt[2][kptno] + (double) l1;
		mag_pw1 = sqrt(pw1_x*pw1_x + pw1_y*pw1_y+ pw1_z*pw1_z);

		for(pw2=0;pw2<npw;pw2++) {
		  h2 = WFK->kg[kptno][pw2][0];
		  k2 = WFK->kg[kptno][pw2][1];
		  l2 = WFK->kg[kptno][pw2][2];

		  pw2_x = WFK->kpt[0][kptno] + (double) h2;
		  pw2_y = WFK->kpt[1][kptno] + (double) k2;
		  pw2_z = WFK->kpt[2][kptno] + (double) l2;
		  mag_pw2 = sqrt(pw2_x*pw2_x + pw2_y*pw2_y+ pw2_z*pw2_z);
		
		  delta_h = h1 - h2;
		  delta_k = k1 - k2;
		  delta_l = l1 - l2;
		  if ((delta_h!=H)||(delta_k!=K)||(delta_l!=L)) {
			continue;
		  }
		  if ((mag_pw1!=mag_pw2)) continue;
		  printf( "\tHKL = %d %d %d\n", delta_h, delta_k, delta_l);
		  printf( "\t\tpw1 = %lf %lf %lf\n", pw1_x, pw1_y, pw1_z);
		  printf( "\t\tpw2 = %lf %lf %lf\n", pw2_x, pw2_y, pw2_z);
		  if (delta_h < 0) delta_h+=ngfftx;
		  if (delta_k < 0) delta_k+=ngffty;
		  if (delta_l < 0) delta_l+=ngfftz;

		  /*loop over bands to determine to asssign coeff and dE*/
		  for(band=0;band<nband;band++) {
			/*Determine Energy range to store contribution*/
			bandE = (((WFK->eigen[kptno][band]-fermi)*HATOEV-bandE_min)*Emesh) + 0.5;
			dE = floor(bandE);
			if (dE >= nEstep) continue; //skip cont if range too high
			
			/*Setting Band Occupation and Kpt Weight*/
			//occ = WFK->occ[kptno][band];
			occ = 1;
			
			/*Calculate potential energy for reflection*/
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
			noncomplex = occ*rflc_mult/UC->bohr_cellV;
			potentialE = gsl_complex_mul_real(complexnums, noncomplex);

			ECON->rflc_local[dE][n] +=  GSL_REAL(potentialE);
			total_local += GSL_REAL(potentialE);
		  } /*END band->nband loop*/
	    } /*END pw2->npw loop*/
	  } /*END pw1->npw loop*/
	  printf( "    kpt# %d\tLocal Potential Energy = %lf\n", n, total_local);
	} /*END: kpt loop*/
  } /*END: npxrd loop*/

  printf( "Total Local Potential Energy = %lf\n", total_local);

  /*free allocated variables*/
  BIN->rec_grid = FreeMemory_threeD_complex(BIN->rec_grid, ngfftx, ngffty);
  BIN->cc_rec_grid = FreeMemory_threeD_complex(BIN->rec_grid, ngfftx, ngffty);

}   //END of Calculate kinetic energy contribution function 
