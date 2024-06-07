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

void FFTon_RealGrid(FILE * flog, BinaryGrid* BIN, NumberGrid* GRD, UnitCell* UC)
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

  printf("\nFFT from Direct Space to Reciprocal Space\n");
  fprintf(flog, "\nFFT from Direct Space to Reciprocal Space\n");
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
  for (jz=0;jz<NGZ;jz++) {
    for(jy=0;jy<NGY;jy++) {
      for(jx=0;jx<NGX;jx++) {
        i_index = jz*ngfftx*ngffty+jy*ngfftx+jx;
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
  fprintf(flog, "000 in Reciprocal space after FFT = %lf\n", c_ooo);

  /*Freeing memory used in FFT*/
  fftw_destroy_plan(grid_plan);
  fftw_free(grid_in); fftw_free(grid_out);
  BIN->real_grid = FreeMemory_threeD_double(BIN->real_grid, NGX, NGY);

} //END of FFTon_RealGrid
