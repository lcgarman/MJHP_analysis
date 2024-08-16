#include <fftw3.h>
#include <complex.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "fftw_functions.h"

void FFTon_RealGrid( BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC)
{
  int ngfftx, ngffty, ngfftz;
  int NGX, NGY, NGZ;
  int ngfft_size;
  int i_index;
  int o_index;
  int jx, jy, jz;
  int h0, k0, l0;
  double c_ooo;
  fftw_complex * grid_in;
  fftw_complex * grid_out;
  fftw_plan grid_plan;
  double re_rec_grid;
  double im_rec_grid;

  printf("\tFFT from Direct Space to Reciprocal Space\n");
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = GRD->NGX;
  NGY = GRD->NGY;
  NGZ = GRD->NGZ;
  /*Allocating Memory for fftw grids*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_BACKWARD, FFTW_MEASURE); 
  BIN->rec_grid = AllocateMemory_threeD_complex(BIN->rec_grid, ngfftx, ngffty, ngfftz);
  BIN->cc_rec_grid = AllocateMemory_threeD_complex(BIN->cc_rec_grid, ngfftx, ngffty, ngfftz);

  /*fill input fft grid with Binary real space grid*/
  for(jx=0;jx<NGX;jx++) {
    for(jy=0;jy<NGY;jy++) {
	  for (jz=0;jz<NGZ;jz++) {
        i_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
	    if ((jx<ngfftx)&&(jy<ngffty)&&(jz<ngfftz)) {
		  grid_in[i_index][REAL] = BIN->real_grid[jx][jy][jz];
		  grid_in[i_index][IMAG] = 0.0;
		} 
	  }
	}
  }
  /*Execute the fft*/
  fftw_execute(grid_plan);

  /*unwrap the fft in reciprocal space*/
  c_ooo = 0.0;
  for(h0=0;h0<ngfftx;h0++) {
    for(k0=0;k0<ngffty;k0++) {
      for(l0=0;l0<ngfftz;l0++) {
        o_index = h0*ngfftz*ngffty+k0*ngfftz+l0;
	    re_rec_grid = grid_out[o_index][REAL]*UC->voxelV;
	    im_rec_grid = grid_out[o_index][IMAG]*UC->voxelV;
		BIN->rec_grid[h0][k0][l0] = gsl_complex_rect(re_rec_grid, im_rec_grid);
		BIN->cc_rec_grid[h0][k0][l0] = gsl_complex_conjugate(BIN->rec_grid[h0][k0][l0]);
 	    if ((h0==0)&&(k0==0)&&(l0==0)) {
		  c_ooo += re_rec_grid;
		}
	  }
	}
  }
  printf( "000 in Reciprocal space after FFT = %lf\n", c_ooo);

  /*Freeing memory used in FFT*/
  fftw_destroy_plan(grid_plan);
  fftw_free(grid_in); fftw_free(grid_out);
  BIN->real_grid = FreeMemory_threeD_double(BIN->real_grid, NGX, NGY);

} //END of FFTon_RealGrid

void FFTon_ReciprocalGrid(BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC, Symmetry* SYM, BinaryGrid* POT)
{
  int ngfftx, ngffty, ngfftz;
  int NGX, NGY, NGZ;
  int ngfft_size;
  int h0, k0, l0;
  int jx, jy, jz;
  int i_index, o_index;
  fftw_complex * grid_in;
  fftw_complex * grid_out;
  fftw_plan grid_plan;
  double re_grid;
  double im_grid;
  double*** real_grid;
  double*** sym_real_grid;
  double coeff_total;
  double F000, f000;
  int symm_error_flag;
  double zf, yf, xf;
  double xf2, yf2, zf2;
  int nsym;
  double ix_d, iy_d, iz_d;
  int stop;
  int ix, iy, iz;
  int j;

  nsym = SYM->nsym;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  NGX = ngfftx+1;
  NGY = ngffty+1;
  NGZ = ngfftz+1;
  coeff_total = 0.0;
  real_grid = NULL;
  sym_real_grid = NULL;

  /*Allocating*/
  real_grid = AllocateMemory_threeD_double(real_grid, ngfftx, ngffty, ngfftz);
  sym_real_grid = AllocateMemory_threeD_double(sym_real_grid, ngfftx, ngffty, ngfftz);
  /*Memory Allocation for fftw3*/
  ngfft_size = ngfftx*ngffty*ngfftz;
  grid_in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ngfft_size);
  grid_plan = fftw_plan_dft_3d(ngfftx, ngffty, ngfftz, grid_in, grid_out, FFTW_FORWARD, FFTW_MEASURE); 
  if ((grid_in==NULL)||(grid_out==NULL)) {
    printf("ERROR: Memory Allocation Failed for wfk to density fftw\n");
    exit(0);
  }
  /*End of memory allocation*/

  printf("\n Performing FFT on Potential to transfrom to real space grid.\n");
  coeff_total = 0.0;
  /*fill fftw 3d grid to perform FFT*/
  for (h0=0;h0<ngfftx;h0++) {
    for (k0=0;k0<ngffty;k0++) {
      for (l0=0;l0<ngfftz;l0++) {
        i_index = h0*ngfftz*ngffty+k0*ngfftz+l0;
		grid_in[i_index][REAL] = GSL_REAL(POT->rec_grid[h0][k0][l0]);
		grid_in[i_index][IMAG] = GSL_IMAG(POT->rec_grid[h0][k0][l0]);
		coeff_total += grid_in[i_index][REAL];
      }
    }
  }
  F000 = GSL_REAL(POT->rec_grid[0][0][0]);
  printf("Total Potential Energy Before FFT: %lf (F000 = %lf)\n", coeff_total, F000);

  /*Execute the fft*/
  fftw_execute(grid_plan);

  coeff_total = 0.0;
  /*unwrap fft in real space*/
  for (jx=0;jx<ngfftx;jx++) {
    for (jy=0;jy<ngffty;jy++) {
      for (jz=0;jz<ngfftz;jz++) {
        o_index = jx*ngfftz*ngffty+jy*ngfftz+jz;
        re_grid = grid_out[o_index][REAL]/UC->bohr_cellV;
        im_grid = grid_out[o_index][IMAG]/UC->bohr_cellV;
        real_grid[jx][jy][jz] = re_grid*UC->voxelV;
        coeff_total += real_grid[jx][jy][jz];
      }
    }
  }
  f000 = real_grid[0][0][0]*ngfft_size;
  printf("Total Potential Energy After FFT: %lf (f000 = %lf)\n", coeff_total, f000);

  fftw_destroy_plan(grid_plan);
  fftw_free(grid_in);
  fftw_free(grid_out);
  POT->rec_grid = FreeMemory_threeD_complex(POT->rec_grid, ngfftx, ngffty);
  
  /*allocate memory for real space potential grid*/
  BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);

  /*symmeterize the potential grid in real space*/
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
  f000 = BIN->real_grid[0][0][0]*ngfft_size;
  printf("Total Symmetrized Potential Energy =  %lf (f000=%lf)\n", coeff_total, f000);

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


