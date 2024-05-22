#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "cell_parameters.h"

void Determine_CellParameters(FILE * flog, UnitCell* UC, NumberGrid* GRD) 
{  
  /*Calculates real and reciprocal space latt params*/
  fprintf(flog, "\n Cell Parameters:\n");
  fprintf(flog, "Real Space Lattice Parameters (bohr): \n");
  fprintf(flog, "a = %lf\t%lf\t%lf\n", UC->bohr_ax, UC->bohr_ay, UC->bohr_az);
  fprintf(flog, "b = %lf\t%lf\t%lf\n", UC->bohr_bx, UC->bohr_by, UC->bohr_bz);
  fprintf(flog, "c = %lf\t%lf\t%lf\n", UC->bohr_cx, UC->bohr_cy, UC->bohr_cz);
  
  /*Determine cell volume in bohr^3 using primitive vectors*/  
  UC->bohr_cellV = UC->bohr_ax*(UC->bohr_by*UC->bohr_cz-UC->bohr_bz*UC->bohr_cy)-\
				   UC->bohr_ay*(UC->bohr_bx*UC->bohr_cz-UC->bohr_bz*UC->bohr_cx)+\
				   UC->bohr_az*(UC->bohr_bx*UC->bohr_cy-UC->bohr_by*UC->bohr_cx);
  //printf("Cell Volume in Bohr^3 = %lf\n", UC->bohr_cellV);

  /*Determine the voxel volume in bohr^3*/
  UC->voxelV = UC->bohr_cellV/(GRD->ngfftx*GRD->ngffty*GRD->ngfftz);
  //printf("Voxel Volume in Bohr^3 = %lf\n", UC->voxelV);

  /*Determine Reciprocal lattice parameters in Bohr*/
  UC->bohr_ax_star = 2*PI*(UC->bohr_by*UC->bohr_cz-UC->bohr_bz*UC->bohr_cy)/UC->bohr_cellV;
  UC->bohr_ay_star = -2*PI*(UC->bohr_bx*UC->bohr_cz-UC->bohr_bz*UC->bohr_cx)/UC->bohr_cellV;
  UC->bohr_az_star = 2*PI*(UC->bohr_bx*UC->bohr_cy-UC->bohr_by*UC->bohr_cx)/UC->bohr_cellV;
  UC->bohr_bx_star = 2*PI*(UC->bohr_cy*UC->bohr_az-UC->bohr_cz*UC->bohr_ay)/UC->bohr_cellV;
  UC->bohr_by_star = -2*PI*(UC->bohr_cx*UC->bohr_az-UC->bohr_cz*UC->bohr_ax)/UC->bohr_cellV;
  UC->bohr_bz_star = 2*PI*(UC->bohr_cx*UC->bohr_ay-UC->bohr_cy*UC->bohr_ax)/UC->bohr_cellV;
  UC->bohr_cx_star = 2*PI*(UC->bohr_ay*UC->bohr_bz-UC->bohr_az*UC->bohr_by)/UC->bohr_cellV;
  UC->bohr_cy_star = -2*PI*(UC->bohr_ax*UC->bohr_bz-UC->bohr_az*UC->bohr_bx)/UC->bohr_cellV;
  UC->bohr_cz_star = 2*PI*(UC->bohr_ax*UC->bohr_by-UC->bohr_ay*UC->bohr_bx)/UC->bohr_cellV;
  //printf("Reciprocal lattice parameters in bohr^-1 =\n");
  //printf("a* = %lf\t%lf\t%lf\n", UC->bohr_ax_star, UC->bohr_ay_star, UC->bohr_az_star);
  //printf("b* = %lf\t%lf\t%lf\n", UC->bohr_bx_star, UC->bohr_by_star, UC->bohr_bz_star);
  //printf("c* = %lf\t%lf\t%lf\n", UC->bohr_cx_star, UC->bohr_cy_star, UC->bohr_cz_star);

  /*Convert everything from units bohr to angstorms*/
  UC->ang_ax=UC->bohr_ax*BOHRTOANG;
  UC->ang_ay=UC->bohr_ay*BOHRTOANG;
  UC->ang_az=UC->bohr_az*BOHRTOANG;
  UC->ang_bx=UC->bohr_bx*BOHRTOANG;
  UC->ang_by=UC->bohr_by*BOHRTOANG;
  UC->ang_bz=UC->bohr_bz*BOHRTOANG;
  UC->ang_cx=UC->bohr_cx*BOHRTOANG;
  UC->ang_cy=UC->bohr_cy*BOHRTOANG;
  UC->ang_cz=UC->bohr_cz*BOHRTOANG;
  fprintf(flog, "Real Space lattice parameters in ang =\n");
  fprintf(flog, "a = %lf\t%lf\t%lf\n", UC->ang_ax, UC->ang_ay, UC->ang_az);
  fprintf(flog, "b = %lf\t%lf\t%lf\n", UC->ang_bx, UC->ang_by, UC->ang_bz);
  fprintf(flog, "c = %lf\t%lf\t%lf\n", UC->ang_cx, UC->ang_cy, UC->ang_cz);

  /*Determine unit cell volume in angstroms*/
  UC->ang_cellV = UC->ang_ax*(UC->ang_by*UC->ang_cz-UC->ang_cy*UC->ang_bz)-\
				  UC->ang_ay*(UC->ang_bx*UC->ang_cz-UC->ang_cx*UC->ang_bz)+\
 				  UC->ang_az*(UC->ang_bx*UC->ang_cy-UC->ang_cx*UC->ang_by);  
  //printf("Cell Volume in ang^3 = %lf\n", UC->ang_cellV);
  /*Determine reciprocal lattice parameters in angstroms*/
  UC->ang_ax_star = (2*PI)*(UC->ang_by*UC->ang_cz-UC->ang_cy*UC->ang_bz)/UC->ang_cellV;
  UC->ang_ay_star = -(2*PI)*(UC->ang_bx*UC->ang_cz-UC->ang_cx*UC->ang_bz)/UC->ang_cellV;
  UC->ang_az_star = (2*PI)*(UC->ang_bx*UC->ang_cy-UC->ang_cx*UC->ang_by)/UC->ang_cellV;
  UC->ang_bx_star = (2*PI)*(UC->ang_cy*UC->ang_az-UC->ang_ay*UC->ang_cz)/UC->ang_cellV;
  UC->ang_by_star = -(2*PI)*(UC->ang_cx*UC->ang_az-UC->ang_ax*UC->ang_cz)/UC->ang_cellV;
  UC->ang_bz_star = (2*PI)*(UC->ang_cx*UC->ang_ay-UC->ang_ax*UC->ang_cy)/UC->ang_cellV;
  UC->ang_cx_star = (2*PI)*(UC->ang_ay*UC->ang_bz-UC->ang_by*UC->ang_az)/UC->ang_cellV;
  UC->ang_cy_star = -(2*PI)*(UC->ang_ax*UC->ang_bz-UC->ang_bx*UC->ang_az)/UC->ang_cellV;
  UC->ang_cz_star = (2*PI)*(UC->ang_ax*UC->ang_by-UC->ang_ay*UC->ang_bx)/UC->ang_cellV;
  fprintf(flog, "Reciprocal lattice parameters in ang^-1 =\n");
  fprintf(flog, "a* = %lf\t%lf\t%lf\n", UC->ang_ax_star, UC->ang_ay_star, UC->ang_az_star);
  fprintf(flog, "b* = %lf\t%lf\t%lf\n", UC->ang_bx_star, UC->ang_by_star, UC->ang_bz_star);
  fprintf(flog, "c* = %lf\t%lf\t%lf\n", UC->ang_cx_star, UC->ang_cy_star, UC->ang_cz_star);
	
}  //END of Determine_CellParameters Function

void Calculate_FermiRadius(UnitCell* UC, FermiSphere* FS) 
{
  double Ghkl_fermi;

  Ghkl_fermi = 2*cbrt(3*PI*PI*FS->vec/UC->ang_cellV);	
  FS->rad_Ef = Ghkl_fermi;
} //END of Calculate_FermiRadius

void Calculate_FermiDegree(UnitCell* UC, FermiSphere* FS) 
{
  double Ghkl_fermi;
  double d_hkl;
  double rad_theta;
  double deg_theta;
  double two_theta;
  double lambda = 0.709319;
  
  Ghkl_fermi = 2*cbrt(3*PI*PI*FS->vec/UC->ang_cellV);	
  d_hkl = (2*PI)/Ghkl_fermi;
  rad_theta = asin((1*lambda)/(2*d_hkl));
  deg_theta = rad_theta*(180/PI);
  two_theta = deg_theta*2;
  FS->two_theta = two_theta;
  //printf("Fermi Radius = %lf (two theta)\n", FS->two_theta);
}

void symmorphic_symmetry(Symmetry* SYM) 
{
  int nsym;
  int sym;
  int nsymor;
  int n;
  int*** symor;
  double tolerance;

  nsym = SYM->nsym; 
  tolerance = 1e-10;
  symor = NULL;

  /*allocate memory for temporary symrel vars*/
  symor = AllocateMemory_threeD_int(symor, 3, 3, nsym);

  nsymor = 0;
  for (sym=0;sym<nsym;sym++) {
	if ((fabs(SYM->tnons[0][sym])>tolerance)||(fabs(SYM->tnons[1][sym])>tolerance)||(fabs(SYM->tnons[2][sym])>tolerance)) continue;
    symor[0][0][nsymor] = SYM->symrel[0][0][sym];
    symor[0][1][nsymor] = SYM->symrel[0][1][sym];
    symor[0][2][nsymor] = SYM->symrel[0][2][sym];
    symor[1][0][nsymor] = SYM->symrel[1][0][sym];
    symor[1][1][nsymor] = SYM->symrel[1][1][sym];
    symor[1][2][nsymor] = SYM->symrel[1][2][sym];
    symor[2][0][nsymor] = SYM->symrel[2][0][sym];
    symor[2][1][nsymor] = SYM->symrel[2][1][sym];
    symor[2][2][nsymor] = SYM->symrel[2][2][sym];
	//printf(" %d\n\t%d %d %d\n\t%d %d %d\n\t%d %d %d\n\n", nsymor, symor[0][0][nsymor], symor[0][1][nsymor], symor[0][2][nsymor], symor[1][0][nsymor], symor[1][1][nsymor], symor[1][2][nsymor], symor[2][0][nsymor], symor[2][1][nsymor], symor[2][2][nsymor]);
    nsymor++;
  }
  
  /*allocate memory for symor is struct SYM*/
  SYM->symor = AllocateMemory_threeD_int(SYM->symor, 3, 3, nsymor);

  SYM->nsymor = nsymor;
  for (n=0;n<nsymor;n++) {
	SYM->symor[0][0][n] = symor[0][0][n];
	SYM->symor[0][1][n] = symor[0][1][n];
	SYM->symor[0][2][n] = symor[0][2][n];
	SYM->symor[1][0][n] = symor[1][0][n];
	SYM->symor[1][1][n] = symor[1][1][n];
	SYM->symor[1][2][n] = symor[1][2][n];
	SYM->symor[2][0][n] = symor[2][0][n];
	SYM->symor[2][1][n] = symor[2][1][n];
	SYM->symor[2][2][n] = symor[2][2][n];
  }
  //printf("# symmetry operations = %d\t # nonsymmorphic = %d\n", SYM->nsym, SYM->nsymor);
  
  /*free temporary symor array*/
  symor = FreeMemory_threeD_int(symor, 3, 3);
 
} /*end of symmorphic symmetry*/
