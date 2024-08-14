#include <fftw3.h>
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
#include "mjhpHKL_density.h"

void mjhpHKL_density(NumberGrid *GRD, BinaryGrid *BIN, Wavefunction *WFK, Symmetry *SYM, UnitCell *UC, VectorIndices *VECT)
{
  int nkpt;
  int nband;
  int npw;
  int kptno;
  int band;
  double bandE;
  double fermi;
  double scanE_start, scanE_stop;
  int pw1;
  double pw1_x, pw1_y, pw1_z;
  double pw2_x, pw2_y, pw2_z;
  double ang_pw1x, ang_pw1y, ang_pw1z;
  double ang_pw2x, ang_pw2y, ang_pw2z;
  double mag_pw1, mag_pw2;
  double minr, maxr;
  int h1, k1, l1;
  int hpos, kpos, lpos;
  int ngfftx, ngffty, ngfftz;
  int NGX, NGY, NGZ;
  int nHKL;
  int MJ_H, MJ_K, MJ_L;
  int HKL_mult;
  double sigma;
  double mag_diff;
  double exponent;
  double broad;
  double wtk;
  gsl_complex c1;
  double coeff_total;
  double*** real_grid;
  double*** sym_real_grid;
  double occ;
  int jx, jy, jz;
  double re_grid;
  double im_grid;
  int i_index;
  int o_index;
  fftw_complex * grid_in;
  fftw_complex * grid_out;
  fftw_plan wfk_den_plan;
  int ngfft_size;
  int symm_error_flag;
  double zf, yf, xf;
  double xf2, yf2, zf2;
  int nsym;
  double ix_d, iy_d, iz_d;
  int stop;
  int ix, iy, iz;
  int j;

  nkpt = WFK->nkpt;
  nband = WFK->nband;
  fermi = WFK->fermi;
  nsym = SYM->nsym;
  minr = VECT->minr;
  maxr = VECT->maxr;
  scanE_start = VECT->scanE_start;
  scanE_stop = VECT->scanE_stop;
  nHKL = VECT->nHKL;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = GRD->ngfftx+1;
  NGY = GRD->ngffty+1;
  NGZ = GRD->ngfftz+1;
  sigma = SIGMA;
  coeff_total = 0.0;
  real_grid = NULL;
  sym_real_grid = NULL;

  /*Allocating for WavefunctionCoefficients*/
  real_grid = AllocateMemory_threeD_double(real_grid, ngfftx, ngffty, ngfftz);
  sym_real_grid = AllocateMemory_threeD_double(sym_real_grid, ngfftx, ngffty, ngfftz);
  /*Memory Allocation for fftw3*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  wfk_den_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_FORWARD, FFTW_MEASURE); 
  if ((grid_in==NULL)||(grid_out==NULL)) {
    printf("ERROR: Memory Allocation Failed for wfk to density fftw\n");
    exit(0);
  }
  /*End of memory allocation*/
  
  printf("\nCalculating Electron Density\n");
  /*Begin Calculating Electron density from HKL*/
  for (kptno=0;kptno<nkpt;kptno++) { 
    printf( "kpt %d \t%lf %lf %lf \tmult=%d\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno], SYM->mult[kptno]);
	npw = WFK->npw[kptno];

	for (band=0;band<nband;band++) {
	  /*Scanning energy within minE eV to maxE eV range*/
	  bandE = (WFK->eigen[kptno][band] - fermi)*HATOEV; 
	  if ((bandE > scanE_stop) || (bandE < scanE_start)) continue; 
      printf("\tband = %d E = %lf\n", band, bandE);
	  
	  HKL_mult = 0;
	  /*zero out grid_in array*/
	  for (h1=0;h1<ngfftx;h1++) {
		for (k1=0;k1<ngffty;k1++) {
		  for (l1=0;l1<ngfftz;l1++) {
		    i_index = h1*ngfftz*ngffty+k1*ngfftz+l1;
		    grid_in[i_index][REAL] = 0.0;
		    grid_in[i_index][IMAG] = 0.0;
		  }
	    }
	  }

	  for (pw1=0;pw1<npw;pw1++) {
		/*Define first planwave coords*/ 
		h1 = WFK->kg[kptno][pw1][0];
		k1 = WFK->kg[kptno][pw1][1];
		l1 = WFK->kg[kptno][pw1][2];
		/*find |k1+Ghkl| */
		pw1_x = WFK->kpt[0][kptno] + (double) h1;
		pw1_y = WFK->kpt[1][kptno] + (double) k1;
		pw1_z = WFK->kpt[2][kptno] + (double) l1;
		ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
		ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
		ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
		mag_pw1 = sqrt(ang_pw1x*ang_pw1x + ang_pw1y*ang_pw1y + ang_pw1z*ang_pw1z);
		/*continue if |pw1|>maxr or |pw1|<minr*/
		if ((mag_pw1>maxr)||(mag_pw1<minr)) continue;
		
		for(j=0;j<nHKL;j++) {
		  MJ_H =  VECT->H_arr[j];
		  MJ_K =  VECT->K_arr[j];
		  MJ_L =  VECT->L_arr[j];
		  /*find pw2 = pw1-HKL*/
		  pw2_x = pw1_x - (double) MJ_H; 
		  pw2_y = pw1_y - (double) MJ_K; 
		  pw2_z = pw1_z - (double) MJ_L; 
		  ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		  ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		  ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
		  mag_pw2 = sqrt(ang_pw2x*ang_pw2x + ang_pw2y*ang_pw2y + ang_pw2z*ang_pw2z);
		  /*if both pw1 and pw2 lie in shell and pw1-pw2=HKL*/
		  if ((mag_pw2>maxr)||(mag_pw2<minr)) continue;
		  HKL_mult++;
		  printf("\tH K L = %d %d %d  mult=%d\n", MJ_H, MJ_K, MJ_L, HKL_mult); 
		  printf("\t\t pw1 = %lf %lf %lf |pw1|=%lf\n", pw1_x, pw1_y, pw1_z, mag_pw1);
		  printf("\t\t pw2 = %lf %lf %lf |pw2|=%lf\n", pw2_x, pw2_y, pw2_z, mag_pw2);
		  /*calculate broadening*/
		  mag_diff = mag_pw1 - mag_pw2;
		  exponent = -(mag_diff*mag_diff)/sigma;
		  broad = exp(exponent);
		  printf( "\t\tbroadening = %lf\n", broad);
	    }
		  
		/*make a hkl grid needs to be all (+) numbers*/
		if (h1<0) hpos = h1 + ngfftx;
		else hpos = h1;
		if (k1<0) kpos = k1 + ngffty;
		else kpos = k1;
		if (l1<0) lpos = l1 + ngfftz;
		else lpos = l1;
		i_index = hpos*ngfftz*ngffty+kpos*ngfftz+lpos;
		  
		/*fill fft grid with wavefunction coefficients*/
		grid_in[i_index][REAL] = WFK->cg[kptno][band][pw1][0];
		grid_in[i_index][IMAG] = -WFK->cg[kptno][band][pw1][1];

	  }
      fftw_execute(wfk_den_plan);

	  /*Setting Band Occupation and Kpt Weight*/
//      occ = WFK->occ[kptno][band];
	  occ = 1;
      wtk = WFK->wtk[kptno];

	  for(jx=0;jx<ngfftx;jx++) {
		for(jy=0;jy<ngffty;jy++) {
		  for(jz=0;jz<ngfftz;jz++) {
			o_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
            re_grid = grid_out[o_index][REAL];
            im_grid = grid_out[o_index][IMAG];
            c1 = gsl_complex_rect(re_grid, im_grid);
            real_grid[jx][jy][jz] += HKL_mult*broad*wtk*occ*gsl_complex_abs2(c1)/UC->bohr_cellV;
            coeff_total += HKL_mult*broad*wtk*occ*gsl_complex_abs2(c1)/UC->bohr_cellV;
		  }  /*END jz->ngfftz loop*/
		}  /*END jy->ngffty loop*/
	  }  /*END jz->ngfftz loop*/
	} /*END band->nband loop*/
    printf("\tTotal Density = %lf\n", coeff_total*UC->voxelV);
  } /*END: kpt loop*/
  printf("\nTotal Mott-Jones Density = %lf\n", coeff_total*UC->voxelV);

  /*free fftw memory*/
  fftw_destroy_plan(wfk_den_plan);
  fftw_free(grid_in);
  fftw_free(grid_out);

  BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);
  /*end of allocation*/

  /*symmeterize the density grid in real space*/
  symm_error_flag = 0;
  for(jx=0;jx<ngfftx;jx++) {
	xf = (double)jx/(double)ngfftx;
	for(jy=0;jy<ngffty;jy++) {
	  yf = (double)jy/(double)ngffty;
	  for(jz=0;jz<ngfftz;jz++) {
		zf = (double)jz/(double)ngfftz;
	    for(j=0;j<nsym;j++) {
	      xf2 = (double)SYM->symrel[0][0][j]*xf + (double)SYM->symrel[0][1][j]*yf + (double)SYM->symrel[0][2][j]*zf + SYM->tnons[0][j];
	      yf2 = (double)SYM->symrel[1][0][j]*xf + (double)SYM->symrel[1][1][j]*yf + (double)SYM->symrel[1][2][j]*zf + SYM->tnons[1][j];
	      zf2 = (double)SYM->symrel[2][0][j]*xf + (double)SYM->symrel[2][1][j]*yf + (double)SYM->symrel[2][2][j]*zf + SYM->tnons[2][j];
	      ix_d = xf2*(double)ngfftx;
	      iy_d = yf2*(double)ngffty;
	      iz_d = zf2*(double)ngfftz;
	      stop = 0;
	      while(stop==0) {
	        if(ix_d < 0.0) ix_d += 1.0*ngfftx;
              else stop = 1;
	      }
	      stop = 0;
	      while(stop==0) {
	        if(iy_d < 0.0) iy_d += 1.0*ngffty;
	        else stop = 1;
	      }
	      stop = 0;
	      while(stop==0) {
	        if(iz_d < 0.0) iz_d += 1.0*ngfftz;
	        else stop = 1;
	      }
	      /*floor function rounds value down to nearest int*/
	      ix = (int)(floor(ix_d+0.5)); 
	      iy = (int)(floor(iy_d+0.5)); 
	      iz = (int)(floor(iz_d+0.5));
	      /*fabs function returns absolute value*/ 
	      if(fabs((double) ix-ix_d) > 0.001) symm_error_flag = 1;
	      if(fabs((double) iy-iy_d) > 0.001) symm_error_flag = 1;
	      if(fabs((double) iz-iz_d) > 0.001) symm_error_flag = 1;
          /*remainder of i/ngfft to ensure i is within range (0,ngfft-1)*/
	      ix = ix % ngfftx; 
	      iy = iy % ngffty;
	      iz = iz % ngfftz;
	      sym_real_grid[ix][iy][iz] += real_grid[jx][jy][jz]/nsym;
	    }
	  }
    }
  }/*end of symm check*/
  if(symm_error_flag==1) {
    printf("Warning: ngfft grid spacing is incompatible with space group symmetry.\n");
  }

  coeff_total = 0.0;
  for(jx=0;jx<ngfftx;jx++) {
    for(jy=0;jy<ngffty;jy++) {
	  for(jz=0;jz<ngfftz;jz++) {
        BIN->real_grid[jx][jy][jz] = sym_real_grid[jx][jy][jz];
        coeff_total += BIN->real_grid[jx][jy][jz];
      }
    }
  }
  printf("Total Normalized MJ Density %lf\n", coeff_total*UC->voxelV);

  /*rewrap grid so jxyz[max] = jxyz[0]*/
  for(jx=0;jx<NGX;jx++) {
    for(jy=0;jy<NGY;jy++) {
	  for(jz=0;jz<NGZ;jz++) {
        if (jx==ngfftx) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][jy][jz];
        if (jy==ngffty) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][0][jz];
        if (jz==ngfftz) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][jy][0];
        if ((jx==ngfftx)&&(jy==ngffty)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][0][jz];
        if ((jx==ngfftx)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][jy][0];
        if ((jy==ngffty)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[jx][0][0];
        if ((jx==ngfftx)&&(jy==ngffty)&&(jz==ngfftz)) BIN->real_grid[jx][jy][jz] = BIN->real_grid[0][0][0];
	  }
	}
  }

  real_grid = FreeMemory_threeD_double(real_grid, ngfftx, ngffty);
  sym_real_grid = FreeMemory_threeD_double(sym_real_grid, ngfftx, ngffty);

} //END of jzf_density function

