#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "read_binary_abinit.h"
#include "allocate_memory.h"

void read_binary_abinit(char filename[200], int option, UnitCell* UC, NumberGrid* GRD, Symmetry* SYM, Wavefunction* WFK, BinaryGrid* BIN, AtomicVariables* ATM) 
{
  /*Reads the Abinit Binary outbut files*/

  FILE * fab; /*Pointer to file fab to read Abinit output file*/
  int k, j; 
  char codvsn [110]; /*abinit version number*/
  char title [132]; /*info about psp*/
  int headform; /*header format version*/
  int fform; 
  int bandtot; /*nkpt*nband*/
  int intxc; /*grid breakdown for exchange-correlation E*/
  int date; /*date of calculation*/
  int ixc; /*index of exchange correlation function*/
  int natom; /*number of atoms*/
  int ngfftx; /*number of grid points for fft*/
  int ngffty;
  int ngfftz;
  int NGX, NGY, NGZ;
  int nkpt; /*number of grid points for kpoint generation*/
  int nspden; /*number of spin-density components*/
  int nspinor; /*number of spinorial components of the wavefunctions*/
  int nsppol; /*number of spin polatization*/
  int nsym; /*number of symmetry operations*/
  int npsp; /*number of pseudopotentials to be read*/
  int ntypat; /*number of types of atoms*/
  int pertcase; /*perturbative DFT?*/
  int usepaw; /*use projector augmented waves method*/
  double ecut; /*energy cutoff*/
  double ecutdg; /*energy cutoff for second grid in PAW*/
  double ecutsm; /*energy cutoff smearing*/
  double ecut_eff; 
  double qptnx; /*q-point re-normalized*/
  double qptny;
  double qptnz;
  double rprimd_ax; /*real space primitive translations*/
  double rprimd_ay;
  double rprimd_az;
  double rprimd_bx;
  double rprimd_by;
  double rprimd_bz;
  double rprimd_cx;
  double rprimd_cy;
  double rprimd_cz;
  double stmbias; /*scanning tunneling microscopy bias voltage*/
  double tphysel; /*temperature (physical) of electrons*/
  double tsmear; /*temp of smearing*/
  double znuclpsp;
  double zionpsp;
  int pspso; /*spin-orbit coupling?*/
  int pspdat; /*revision date*/
  int pspcod; 
  int pspxc; /*XC fxnal for psp*/
  int lmn_size; /*spherical harmonics?*/
  int usewvl; /*use wavelet basis set*/
  int istwfkv;
  int nbandv;
  int npwv;
  int so_pspv; /*spint-orbit treatment for each psp should be by NATOM*/
  int symafmv; /*symmetries, anti-ferromagnetic characteristics by MAX_SYM*/
  int symrel_mx, symrel_my, symrel_mz; /*symmetry in real space*/
  int symrel_nx, symrel_ny, symrel_nz;
  int symrel_px, symrel_py, symrel_pz;
  int typatv; /*number of types of atoms*/
  int type;
  double kptx, kpty, kptz;
  double occopt; /*occupation numbers*/
  double tnons_x, tnons_y, tnons_z; /*translation non-symmorphic vect (MAXSYM)*/
  double znucltypatv; /*atomic number for nucs by NATOMS*/
  double* znucltypat; /*atomic number for nucs by NATOMS*/
  double wtkv;
  double residm; /*residual density*/
  double x,y,z;
  double etotal; /*total E*/
  double fermi; /*fermi E as determined by the calc*/
  
  /*Variables for binary real space grid*/
  int jx, jy, jz;
  double bin_grid;

  /*Variable for Symmetry Allocation*/
  double multiplicity;
  int mult_tot;

  /*Variables for Reading wavefunction*/
  int kptno;
  int npw;
  int nband;
  int pw;
  int kx, ky, kz;
  int band;
  double eigenv;  
  double occv;
  double cgv;

  printf("\nReading %s\n", filename);
  /*Open Abinit _o_* output file*/
  fab=fopen(filename,"rb+"); /*open the binary file in read and write mode*/
  if(fab==NULL) {
    printf("ERROR: %s not found. \n", filename);
    exit(0);
  } 

  /*Initializing variables and arrays*/
  j=0; 
  mult_tot = 0;
  znucltypat = NULL;

  //BEGIN READING HEADER
  //First Block
  fread(&j, sizeof(int), 1, fab);
  j=fread(codvsn, sizeof(char), 6, fab);
  fread(&headform, sizeof(int), 1, fab);
  fread(&fform, sizeof(int), 1, fab);
  fread(&j, sizeof(int), 1, fab);

  //Second Block
  fread(&j, sizeof(int), 1, fab);
  fread(&bandtot, sizeof(int), 1, fab);
  fread(&date, sizeof(int), 1, fab);
  fread(&intxc, sizeof(int), 1, fab);
  fread(&ixc, sizeof(int), 1, fab);
  fread(&natom, sizeof(int), 1, fab);
    ATM->natom = natom;
  fread(&ngfftx, sizeof(int), 1, fab);
  fread(&ngffty, sizeof(int), 1, fab);
  fread(&ngfftz, sizeof(int), 1, fab);
  /*Storing ngfft into grid info*/
	GRD->ngfftx=ngfftx;
	GRD->ngffty=ngffty;
	GRD->ngfftz=ngfftz;
	NGX=ngfftx+1;
	NGY=ngffty+1;
	NGZ=ngfftz+1;
	GRD->NGX=NGX;
	GRD->NGY=NGY;
	GRD->NGZ=NGZ;

  fread(&nkpt, sizeof(int), 1, fab);
    WFK->nkpt = nkpt;
  fread(&nspden, sizeof(int), 1, fab);
  fread(&nspinor, sizeof(int), 1, fab);
  fread(&nsppol, sizeof(int), 1, fab);
  fread(&nsym, sizeof(int), 1, fab);
    SYM->nsym = nsym;
  fread(&npsp, sizeof(int), 1, fab);
  fread(&ntypat, sizeof(int), 1, fab);
    ATM->ntypat = ntypat;
  fread(&occopt, sizeof(int), 1, fab);
  fread(&pertcase, sizeof(int), 1, fab);
  fread(&usepaw, sizeof(int), 1, fab);
  fread(&ecut, sizeof(double), 1, fab);
  fread(&ecutdg, sizeof(double), 1, fab);
  fread(&ecutsm, sizeof(double), 1, fab);
  fread(&ecut_eff, sizeof(double), 1, fab);
  fread(&qptnx, sizeof(double), 1, fab);
  fread(&qptny, sizeof(double), 1, fab);
  fread(&qptnz, sizeof(double), 1, fab);
  fread(&rprimd_ax, sizeof(double), 1, fab);
  fread(&rprimd_ay, sizeof(double), 1, fab);
  fread(&rprimd_az, sizeof(double), 1, fab);
  fread(&rprimd_bx, sizeof(double), 1, fab);
  fread(&rprimd_by, sizeof(double), 1, fab);
  fread(&rprimd_bz, sizeof(double), 1, fab);
  fread(&rprimd_cx, sizeof(double), 1, fab);
  fread(&rprimd_cy, sizeof(double), 1, fab);
  fread(&rprimd_cz, sizeof(double), 1, fab);
  /*Store primitive vectors in CellInfor structure*/
	UC->bohr_ax = rprimd_ax;
	UC->bohr_bx = rprimd_bx;
	UC->bohr_cx = rprimd_cx;
	UC->bohr_ay = rprimd_ay;
	UC->bohr_by = rprimd_by;
	UC->bohr_cy = rprimd_cy;
	UC->bohr_az = rprimd_az;
	UC->bohr_bz = rprimd_bz;
	UC->bohr_cz = rprimd_cz;

  fread(&stmbias, sizeof(double), 1, fab);
  fread(&tphysel, sizeof(double), 1, fab);
  fread(&tsmear, sizeof(double), 1, fab);
  fread(&usewvl, sizeof(int), 1, fab);
  fread(&j, sizeof(int), 1, fab);
  fread(&j, sizeof(int), 1, fab);

  /*allocate variables for next section of readin*/
  printf( "\tAllocating Memory for Header Variables:\n");
  WFK->npw = AllocateMemory_oneD_int(WFK->npw, nkpt);
  printf( "\t\tSuccess for npw (x nkpt).\n");
  SYM->symrel = AllocateMemory_threeD_int(SYM->symrel, 3, 3, nsym);
  printf( "\t\tSuccess for symrel (x 3 x 3 x nsym).\n");
  ATM->typat = AllocateMemory_oneD_int(ATM->typat, natom);
  printf( "\t\tSuccess for typat (x natom).\n");
  WFK->kpt = AllocateMemory_twoD_double(WFK->kpt, 3, nkpt);
  printf( "\t\tSuccess for kpt (x 3 x nkpt).\n");
  SYM->tnons = AllocateMemory_twoD_double(SYM->tnons, 3, nsym);
  printf( "\t\tSuccess for tnons (x 3 x nsym).\n");
  znucltypat = AllocateMemory_oneD_double(znucltypat, ntypat);
  printf( "\t\tSuccess for znucltypat (x ntypat).\n");
  ATM->atomicno = AllocateMemory_oneD_int(ATM->atomicno, natom);
  printf( "\t\tSuccess for atomicno (x natom).\n");
  WFK->wtk = AllocateMemory_oneD_double(WFK->wtk, nkpt);
  printf( "\t\tSuccess for wtk (x nkpt).\n");
  SYM->mult = AllocateMemory_oneD_int(SYM->mult, nkpt);
  printf( "\t\tSuccess for mult (x nkpt).\n");
  ATM->xred = AllocateMemory_twoD_double(ATM->xred, 3, natom);
  printf( "\t\tSuccess for xred (x 3 x natom).\n");
  /*end of allocation*/

  //Third Block
  for(j=0;j<(nkpt);j++) {
    fread(&istwfkv, sizeof(int), 1, fab);
  }
  for(j=0;j<(nkpt*nsppol);j++) {
    fread(&nbandv, sizeof(int), 1, fab);
    WFK->nband = nbandv;
  }

int max_npw;
max_npw = 0;
  for(j=0;j<(nkpt);j++) {
    fread(&npwv, sizeof(int), 1, fab);
if(npwv>max_npw) max_npw = npwv;
    WFK->npw[j] = npwv;
  }

  for(j=0;j<(npsp);j++) {
    fread(&so_pspv, sizeof(int), 1, fab);
	//so_psp[j] = so_pspv;
  }
  for(j=0;j<(nsym);j++) {
    fread(&symafmv, sizeof(int), 1, fab);
    //symafm[j] = symafmv;
  }
  for(j=0;j<(nsym);j++) {
    fread(&symrel_mx, sizeof(int), 1, fab);
	SYM->symrel[0][0][j] = symrel_mx;
    fread(&symrel_nx, sizeof(int), 1, fab);
	SYM->symrel[1][0][j] = symrel_nx;
    fread(&symrel_px, sizeof(int), 1, fab);
	SYM->symrel[2][0][j] = symrel_px;
    fread(&symrel_my, sizeof(int), 1, fab);
	SYM->symrel[0][1][j] = symrel_my;
    fread(&symrel_ny, sizeof(int), 1, fab);
	SYM->symrel[1][1][j] = symrel_ny;
    fread(&symrel_py, sizeof(int), 1, fab);
	SYM->symrel[2][1][j] = symrel_py;
    fread(&symrel_mz, sizeof(int), 1, fab);
	SYM->symrel[0][2][j] = symrel_mz;
    fread(&symrel_nz, sizeof(int), 1, fab);
	SYM->symrel[1][2][j] = symrel_nz;
    fread(&symrel_pz, sizeof(int), 1, fab);
	SYM->symrel[2][2][j] = symrel_pz;
//printf(" %d\n\t%d %d %d\n\t%d %d %d\n\t%d %d %d\n\n", j, SYM->symrel[0][0][j], SYM->symrel[0][1][j], SYM->symrel[0][2][j], SYM->symrel[1][0][j], SYM->symrel[1][1][j], SYM->symrel[1][2][j], SYM->symrel[2][0][j], SYM->symrel[2][1][j], SYM->symrel[2][2][j]);
  }
  for(j=0;j<(natom);j++) {
    fread(&typatv, sizeof(int), 1, fab);
      ATM->typat[j] = typatv;
  }

  for(k=0;k<(nkpt);k++) {
	  //reduced coordinates of kpts
    fread(&kptx, sizeof(double), 1, fab);
    fread(&kpty, sizeof(double), 1, fab);
    fread(&kptz, sizeof(double), 1, fab);
    WFK->kpt[0][k] = kptx;
    WFK->kpt[1][k] = kpty;
    WFK->kpt[2][k] = kptz;
  }

  for(j=0;j<(bandtot);j++) {
    fread(&occv, sizeof(double), 1, fab);
  }

  for(j=0;j<(nsym);j++) {
    fread(&tnons_x, sizeof(double), 1, fab);
	SYM->tnons[0][j]=tnons_x;
    fread(&tnons_y, sizeof(double), 1, fab);
	SYM->tnons[1][j]=tnons_y;
    fread(&tnons_z, sizeof(double), 1, fab);
	SYM->tnons[2][j]=tnons_z;
    //printf("%d\n\t%lf %lf %lf\n", j, tnons_x, tnons_y, tnons_z);
  }
  for(j=0;j<(ntypat);j++) {
    fread(&znucltypatv, sizeof(double), 1, fab);
    znucltypat[j] = znucltypatv;
  }
  /*store atomic nuber by atoms*/
  for (j=0;j<natom;j++) {
    type = ATM->typat[j] - 1;
    ATM->atomicno[j] = floor(znucltypat[type]);
  }
  /*Read in the kpt-weight find multiplicity of kpts*/
  for(j=0;j<(nkpt);j++) {
    fread(&wtkv, sizeof(double), 1, fab);
    WFK->wtk[j] = wtkv;
	multiplicity = (WFK->wtk[j]/WFK->wtk[0]);
    SYM->mult[j] = ceil(multiplicity);
	mult_tot += SYM->mult[j];
	//printf("kpt: %lf %lf %lf\tmult = %d\t%d\n", WFK->kpt[0][j], WFK->kpt[1][j], WFK->kpt[2][j], SYM->mult[j], mult_tot);
  }
  SYM->mult_tot = mult_tot;
  printf("Total nkpts = %d\n",  mult_tot);
  fread(&j, sizeof(int), 1, fab);

  for(k=0; k<(npsp);k++) {
    fread(&j, sizeof(int), 1, fab); 
    fread(&title, sizeof(char), 132, fab);
    fread(&znuclpsp, sizeof(double), 1, fab);
    fread(&zionpsp, sizeof(double), 1, fab);
    fread(&pspso, sizeof(int), 1, fab);
    fread(&pspdat, sizeof(int), 1, fab);
    fread(&pspcod, sizeof(int), 1, fab);
    fread(&pspxc, sizeof(int), 1, fab);
    fread(&lmn_size, sizeof(int), 1, fab);
    fread(&j, sizeof(int), 1, fab);
  }

  //Fifth Block
  if(usepaw==0) {
    fread(&j, sizeof(int), 1, fab);
    fread(&residm, sizeof(double), 1, fab);
     /*Xred data*/
    for(k=0;k<natom;k++) {
      fread(&x, sizeof(double), 1, fab);
      fread(&y, sizeof(double), 1, fab);
      fread(&z, sizeof(double), 1, fab);
      ATM->xred[0][k] = x;
      ATM->xred[1][k] = y;
      ATM->xred[2][k] = z;
        /*pointers to the X/Y/Z cart info in struct file */
    }
    fread(&etotal, sizeof(double), 1, fab);
    fread(&fermi, sizeof(double), 1, fab);
    WFK->fermi = fermi;
    fread(&j, sizeof(int), 1, fab);  
  }
  //Sixth Block
  else {
    printf("No code written for usepaw=1. check abinit header for additional info");
  }
  //END of header

  /*Free local variables*/
  znucltypat = FreeMemory_oneD_double(znucltypat);
  
  /*Option=2 => Read only header then exit*/  
  if (option == 2) {
    return;
  }

  /*Option=1 => Read in Binary POT/DEN info*/
  else if (option == 1) {
	printf( " \tAllocating Memory for Binary Grid\n");
    BIN->real_grid = AllocateMemory_threeD_double(BIN->real_grid, NGX, NGY, NGZ);
    printf( "\t\tSuccess for real potential grid (NGX x NGY x NGZ).\n");
    /*Begin reading in binary grid in real space*/
    fread(&j, sizeof(int), 1, fab);
    for(jz=0;jz<NGZ;jz++) {
      for(jy=0;jy<NGY;jy++) {
        for(jx=0;jx<NGX;jx++) {
          if((jx<NGX-1)&&(jy<NGY-1)&&(jz<NGZ-1)) {
            fread(&bin_grid,sizeof(double),1,fab);
            BIN->real_grid[jx][jy][jz]=bin_grid;
          }
        }
      }
    }
    fread(&j, sizeof(int), 1, fab);
    /*END of reading in information from binary file*/
    fclose(fab);

    /*Now unwrap grid so point (jz=ngfftx)==(jz=0)*/
    for(jz=0;jz<NGZ;jz++) {
      for(jy=0;jy<NGY;jy++) {
        for(jx=0;jx<NGX;jx++) {
          if(jx==NGX-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][jy][jz];
          if(jy==NGY-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][0][jz];
          if(jz==NGZ-1) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][jy][0];
          if((jx==NGX-1)&&(jy==NGY-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][0][jz];
          if((jx==NGX-1)&&(jz==NGZ-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][jy][0];
          if((jy==NGY-1)&&(jz==NGZ-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[jx][0][0];
          if((jy==NGY-1)&&(jz==NGZ-1)&&(jx==NGX-1)) BIN->real_grid[jx][jy][jz]=BIN->real_grid[0][0][0];
		}
	  }
	}
	/*End of reading in and manipulating Binary Real space real_grid*/
    printf( "Finished reading %s.\n", filename);
  } //END of if option==1

  /*if option!=1 => Read in WFK file*/
  else {

	/*Allocate Memory for wavefunction variables*/
    printf( "\tAllocating Memory for Wavefunction Variables:\n");
    WFK->eigen = AllocateMemory_twoD_double(WFK->eigen, nkpt, WFK->nband);
    printf( "\t\tSuccess for eigen (nkpt x nband).\n");
    WFK->occ = AllocateMemory_twoD_double(WFK->occ, nkpt, WFK->nband);
    printf( "\t\tSuccess for occ (nkpt, nband).\n");
  fflush(stdout);
    AllocateMemory_Wavefunctions(WFK);
    printf( "\t\tSuccess for cg and kg.\n");
    
	//Being reading WFK info
    for (k=0;k<nsppol;k++) {
      for (kptno=0;kptno<nkpt;kptno++) {
		fflush(stdout);
        fread(&j, sizeof(int), 1, fab);
        fread(&npw, sizeof(int), 1, fab);
        fread(&nspinor, sizeof(int), 1, fab);
        ATM->nspinor = nspinor;
        fread(&nband, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        for (pw=0;pw<npw;pw++) {
          fread(&kx, sizeof(int), 1, fab);
          fread(&ky, sizeof(int), 1, fab);
          fread(&kz, sizeof(int), 1, fab);
	  	  //kxyz = planewave reduced coords
          WFK->kg[kptno][pw][0]=kx;
          WFK->kg[kptno][pw][1]=ky;
          WFK->kg[kptno][pw][2]=kz;
        }
		fread(&j, sizeof(int), 1, fab);
        fread(&j, sizeof(int), 1, fab);
        for (band=0;band<nband;band++) {
          fread(&eigenv, sizeof(double), 1, fab);
	  	  WFK->eigen[kptno][band]=eigenv;
		}
        for (band=0;band<nband;band++) {
          fread(&occv, sizeof(double), 1, fab);
		  WFK->occ[kptno][band]=occv;
        }
		fread(&j, sizeof(int), 1, fab);

		for(band=0;band<nband;band++) {
		  fread(&j, sizeof(int), 1, fab);
          for(pw=0;pw<npw;pw++) {
            fread(&cgv, sizeof(double), 1, fab);
            WFK->cg[kptno][band][pw][0]=cgv;
            fread(&cgv, sizeof(double), 1, fab);
            WFK->cg[kptno][band][pw][1]=cgv;
		  }
		  fread(&j, sizeof(int), 1, fab);

		} //END nband loop 
	  } //END nkpt loop
	}  //END nspoll loop

	/*END of reading in WFK file*/
  fclose(fab);
  printf( "Finished reading %s.\n", filename);
  } //END of else stmt for WFK files

} //END Read_About function

