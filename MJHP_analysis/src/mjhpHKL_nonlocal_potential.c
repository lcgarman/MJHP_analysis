#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_complex_math.h>
#include <gsl/gsl_complex.h>
#include <gsl/gsl_sf_legendre.h>
#include <gsl/gsl_sf_gamma.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "psp_functions.h"
#include "mjhpHKL_nonlocal_potential.h"

void mjhpHKL_nonlocal_potential(FILE * flog, VectorIndices *VECT, NumberGrid* GRD, EnergyStep* ESTP, Wavefunction* WFK, Pseudopotential* PSP, UnitCell* UC, AtomicVariables* ATM, EnergyContribution *ECON)
{ 
  int j;
  int nkpt, kptno;
  int nband, band;
  int npw, pw;
  int PW_max;
  int dE;
  double total_nonlocal;
  int nEstep;
  double bandE_min;
  double fermi;
  double bandE;
  int Emesh;
  double bohr_cellV;
  int ntypat;
  int natom;
  double occ;
  int ngfftx, ngffty, ngfftz;
  int h1, k1, l1;
  int h2, k2, l2;
  int delta_h, delta_k, delta_l;
  int pw1, pw2;
  double* gx;
  double* gy;
  double* gz;
  double* ga;
  double* gb;
  double* gc;
  double* g;
  double* costheta;
  double* phi;
  double* theta;
  double**** p;
  double*** Plm;
  int l, m, i;
  int* l_max_type;
  int** i_max_type;
  int atom_type;
  double**** ph;
  double Dga, Dgb, Dgc;
  double Dphi;
  double theta1, theta2, theta3;
  int table_entry, table_entry1, table_entry2, table_entry3;
  double cos_mdphi_times_2[4];
  double cos_table[10000000];
  double sin_table[10000000];
  double theta_step = 0.000001;
  double* c1c2_RE;
  double* c1c2_IM;
  double php;
  double** vg1g2_real_l;
  double atom_phase_rad;
  double atomphase_RE, atomphase_IM;
  int atomno;
  double gkk_real;
  double c1c2atomphases_real; 
  double pw1_x, pw1_y, pw1_z;
  double pw2_x, pw2_y, pw2_z;
  double ang_pw1x, ang_pw1y, ang_pw1z;
  double ang_pw2x, ang_pw2y, ang_pw2z;
  double mag_pw1, mag_pw2;
  double mag_diff;
  double exponent;
  double broad;
  double minr, maxr;
  int nHKL;
  int hkl_match;
  int H_match, K_match, L_match;
  double wtk;
  double sigma;

  theta1=0;
  for(j=0;j<9999999;j++) {
  cos_table[j]=cos(TWO_PI*theta1);
	sin_table[j]=sin(TWO_PI*theta1);
	theta1+=theta_step;
  }
  
  nkpt = WFK->nkpt;
  nband = WFK->nband;
  ngfftx = GRD->ngfftx;
  ngffty = GRD->ngffty;
  ngfftz = GRD->ngfftz;
  nEstep = ESTP->nEstep;
  bandE_min = ESTP->bandE_min;
  fermi  = WFK->fermi;
  bohr_cellV = UC->bohr_cellV;
  ntypat = ATM->ntypat;
  natom = ATM->natom;
  minr = VECT->minr;
  maxr = VECT->maxr;
  nHKL = VECT->nHKL;
  sigma = SIGMA;
  Emesh = EMESH;
  total_nonlocal = 0.0;

  PW_max = 0;
  for (kptno=0;kptno<nkpt;kptno++) {
    if (WFK->npw[kptno]>PW_max) PW_max = WFK->npw[kptno];
  }

  /*Initialize variables to be allocated*/
  l_max_type = NULL;
  i_max_type = NULL;
  ph = NULL;
  c1c2_RE = NULL;
  c1c2_IM = NULL;
  vg1g2_real_l = NULL;
  gx = NULL;
  gy = NULL;
  gz = NULL;
  ga = NULL;
  gb = NULL;
  gc = NULL;
  g = NULL;
  costheta = NULL;
  phi = NULL;
  theta = NULL;
  p = NULL;
  Plm = NULL;
  /*Allocate memory for large variables in this function*/
  ECON->nonlocal = AllocateMemory_oneD_double(ECON->nonlocal, nEstep);
  l_max_type = AllocateMemory_oneD_int(l_max_type, ntypat);
  i_max_type = AllocateMemory_twoD_int(i_max_type, ntypat, 5);
  ph = AllocateMemory_fourD_double(ph, 4, 4, 5, ntypat);
  c1c2_RE = AllocateMemory_oneD_double(c1c2_RE, nband);
  c1c2_IM = AllocateMemory_oneD_double(c1c2_IM, nband);
  vg1g2_real_l = AllocateMemory_twoD_double(vg1g2_real_l, ntypat, 5);
  gx = AllocateMemory_oneD_double(gx, PW_max);  
  gy = AllocateMemory_oneD_double(gy, PW_max);  
  gz = AllocateMemory_oneD_double(gz, PW_max);  
  ga = AllocateMemory_oneD_double(ga, PW_max);  
  gb = AllocateMemory_oneD_double(gb, PW_max);  
  gc = AllocateMemory_oneD_double(gc, PW_max);  
  g = AllocateMemory_oneD_double(g, PW_max);  
  costheta = AllocateMemory_oneD_double(costheta, PW_max);
  phi = AllocateMemory_oneD_double(phi, PW_max);
  theta = AllocateMemory_oneD_double(theta, PW_max);
  p = AllocateMemory_fourD_double(p, PW_max, 4, 4, ntypat);
  Plm = AllocateMemory_threeD_double(Plm, PW_max, 4, 4);

  /*END of varibale allocation*/

  if (ATM->nspinor!=1) {
    printf("ERROR: No code written for nspinor!=0\n");
    exit(0);
  }
  
  printf("\nCalculating Nonlocal Potential Energy \n");
  /*begin calculating nonlocal energy for HKL*/
  for (kptno=0;kptno<nkpt;kptno++){ 
    fprintf(flog, "kpt %d \t%lf %lf %lf\n", kptno, WFK->kpt[0][kptno], WFK->kpt[1][kptno], WFK->kpt[2][kptno]);
	npw = WFK->npw[kptno];

	/*Find |G+k| needed to calculate kinetic Energy*/
	for(pw=0;pw<npw;pw++) {
	  //gabc = G + k; rec lattice vector (G) plus wavevector(k)
	  ga[pw] = WFK->kpt[0][kptno] + WFK->kg[kptno][pw][0]; 
	  gb[pw] = WFK->kpt[1][kptno] + WFK->kg[kptno][pw][1]; 
	  gc[pw] = WFK->kpt[2][kptno] + WFK->kg[kptno][pw][2]; 
	  //converting to primitive cartensian coordinates 
	  gx[pw] = ga[pw]*UC->bohr_ax_star+gb[pw]*UC->bohr_bx_star+gc[pw]*UC->bohr_cx_star;
	  gy[pw] = ga[pw]*UC->bohr_ay_star+gb[pw]*UC->bohr_by_star+gc[pw]*UC->bohr_cy_star;
	  gz[pw] = ga[pw]*UC->bohr_az_star+gb[pw]*UC->bohr_bz_star+gc[pw]*UC->bohr_cz_star;
	  xyz2sph(gx[pw],gy[pw],gz[pw],&g[pw], &theta[pw], &phi[pw]);
	  costheta[pw] = cos(theta[pw]);
	  /*Populate double Plm array*/
	  for(l=0;l<4;l++) {
		for(m=0;m<l+1;m++) {
		  /*Calculating Plm with legendre polynomial*/
		  Plm[pw][l][m] = gsl_sf_legendre_sphPlm(l, m, costheta[pw]);
		}
	  }
	  for(l=0;l<4;l++) {
		for(atom_type=0;atom_type<ntypat;atom_type++) {
		  l_max_type[atom_type]=0;
		  i_max_type[atom_type][l]=0;
		}
	  }
	  for(l=0;l<4;l++) {
		for(i=0;i<3;i++) {
		  for(atom_type=0;atom_type<ntypat;atom_type++) {
			/*Calculating projectors in rec space*/
			p[pw][l][i][atom_type]=projp(l,i,g[pw],PSP,atom_type+1,bohr_cellV);
			/*if p doesnot equal 0, store 1 and i*/
			if (p[pw][l][i][atom_type]!=0.0) {
			  l_max_type[atom_type]=l+1;
			  i_max_type[atom_type][l]=i+1;
			}
		  }
		}
	  }
	}  //END of pw->npw loop*/

	/*Defining pw range to loop over*/
	for(pw1=0;pw1<npw;pw1++) {
	  h1 = WFK->kg[kptno][pw1][0];
	  k1 = WFK->kg[kptno][pw1][1];
	  l1 = WFK->kg[kptno][pw1][2];
	  /*Find magnitude of pw1*/ 
	  pw1_x = WFK->kpt[0][kptno] + (double) h1;
	  pw1_y = WFK->kpt[1][kptno] + (double) k1;
	  pw1_z = WFK->kpt[2][kptno] + (double) l1;
	  ang_pw1x = pw1_x*UC->ang_ax_star + pw1_y*UC->ang_bx_star + pw1_z*UC->ang_cx_star;
	  ang_pw1y = pw1_x*UC->ang_ay_star + pw1_y*UC->ang_by_star + pw1_z*UC->ang_cy_star;
	  ang_pw1z = pw1_x*UC->ang_az_star + pw1_y*UC->ang_bz_star + pw1_z*UC->ang_cz_star;
	  mag_pw1 = sqrt(ang_pw1x*ang_pw1x+ang_pw1y*ang_pw1y+ang_pw1z*ang_pw1z);
      /*continue if |pw1|>maxr or |pw1|<minr*/
      if ((mag_pw1>maxr)||(mag_pw1<minr)) continue;
	  
	  for(atom_type=0;atom_type<ntypat;atom_type++) {
		for(l=0;l<l_max_type[atom_type];l++) {
		  ph[0][0][l][atom_type] = p[pw1][l][0][atom_type]*PSP->h[0][0][l][atom_type];
		  ph[0][1][l][atom_type] = p[pw1][l][0][atom_type]*PSP->h[0][1][l][atom_type];
		  ph[0][2][l][atom_type] =  p[pw1][l][0][atom_type]*PSP->h[0][2][l][atom_type];
		  ph[1][0][l][atom_type] =  p[pw1][l][1][atom_type]*PSP->h[1][0][l][atom_type];
		  ph[1][1][l][atom_type] =  p[pw1][l][1][atom_type]*PSP->h[1][1][l][atom_type];
		  ph[1][2][l][atom_type] =  p[pw1][l][1][atom_type]*PSP->h[1][2][l][atom_type];
		  ph[2][0][l][atom_type] =  p[pw1][l][2][atom_type]*PSP->h[2][0][l][atom_type];
		  ph[2][1][l][atom_type] =  p[pw1][l][2][atom_type]*PSP->h[2][1][l][atom_type];
		  ph[2][2][l][atom_type] =  p[pw1][l][2][atom_type]*PSP->h[2][2][l][atom_type];
		}
	  }
	  for(pw2=0;pw2<npw;pw2++) {
		/*Defining k2 from pw2*/
		h2 = WFK->kg[kptno][pw2][0];
		k2 = WFK->kg[kptno][pw2][1];
		l2 = WFK->kg[kptno][pw2][2];
	    /*Find magnitude of pw2*/ 
		pw2_x = WFK->kpt[0][kptno] + (double) h2;
		pw2_y = WFK->kpt[1][kptno] + (double) k2;
		pw2_z = WFK->kpt[2][kptno] + (double) l2;
		//mag_pw2 = sqrt(pw2_x*pw2_x + pw2_y*pw2_y+ pw2_z*pw2_z);
		ang_pw2x = pw2_x*UC->ang_ax_star + pw2_y*UC->ang_bx_star + pw2_z*UC->ang_cx_star;
		ang_pw2y = pw2_x*UC->ang_ay_star + pw2_y*UC->ang_by_star + pw2_z*UC->ang_cy_star;
		ang_pw2z = pw2_x*UC->ang_az_star + pw2_y*UC->ang_bz_star + pw2_z*UC->ang_cz_star;
		mag_pw2 = sqrt(ang_pw2x*ang_pw2x+ang_pw2y*ang_pw2y+ang_pw2z*ang_pw2z);
        /*continue if |pw2|>maxr or |pw2|<minr*/
        if ((mag_pw2>maxr)||(mag_pw2<minr)) continue;
		
        /*find difference between pw1 and pw2*/
		delta_h = h1 - h2;
		delta_k = k1 - k2;
		delta_l = l1 - l2;
        hkl_match = 0;
        for(j=0;j<nHKL;j++) {
		  if ((delta_h==VECT->H_arr[j])&&(delta_k==VECT->K_arr[j])&&(delta_l==VECT->L_arr[j])) {
            hkl_match = 1;
            H_match = VECT->H_arr[j];
            K_match = VECT->K_arr[j];
            L_match = VECT->L_arr[j];
            break;
		  }
        }
        /*if no match is found then go to next pw2*/
        if (hkl_match == 0) continue;
        fprintf(flog, "\t H K L = %d %d %d\n", H_match, K_match, L_match);
        fprintf(flog, "\t\t pw1 = %lf %lf %lf |pw1|=%lf\n", pw1_x, pw1_y, pw1_z, mag_pw1);
        fprintf(flog, "\t\t pw2 = %lf %lf %lf |pw2|=%lf\n", pw2_x, pw2_y, pw2_z, mag_pw2);
        /*make delta_h positive*/
		if (delta_h < 0) delta_h+=ngfftx;
		if (delta_k < 0) delta_k+=ngffty;
		if (delta_l < 0) delta_l+=ngfftz;

        /*calculate broadening*/
        mag_diff = mag_pw1 - mag_pw2;
        exponent = -(mag_diff*mag_diff)/sigma;
        broad = exp(exponent);
        fprintf(flog, "\t\tbroadening = %lf\n", broad);

		/*Finding delta(reduced plane wave coords in cart*/
		Dga=ga[pw1]-ga[pw2];
		Dgb=gb[pw1]-gb[pw2];
		Dgc=gc[pw1]-gc[pw2];
		Dphi=(phi[pw1]-phi[pw2])*ONE_TWO_PI;
		/*Define theta from phi*/
		theta1 = Dphi;
		theta2 = 2*Dphi;
		theta3 = 3*Dphi;
		table_entry1 = theta1;
		theta1 -= table_entry1;
		table_entry2 = theta2;
		theta2 -= table_entry2;
		table_entry3 = theta3;
		theta3 -= table_entry3;
		if(theta1 < 0.0) {
		  theta1+=1.00;
		  theta2+=1.00;
		  theta3+=1.00;
		}
		/*Assign the (int)table entries based off of (doub) theta vals*/
		table_entry1 = theta1/theta_step + 0.5;
		table_entry2 = theta2/theta_step + 0.5;
		table_entry3 = theta3/theta_step + 0.5;
		cos_mdphi_times_2[1]=2.0*cos_table[table_entry1];
		cos_mdphi_times_2[2]=2.0*cos_table[table_entry2];
		cos_mdphi_times_2[3]=2.0*cos_table[table_entry3];
		/*No code written for nspinor!=0*/
		for(band=0;band<nband;band++) {
		  c1c2_RE[band] = WFK->cg[kptno][band][pw1][0]*WFK->cg[kptno][band][pw2][0] + WFK->cg[kptno][band][pw1][1]*WFK->cg[kptno][band][pw2][1];
		  c1c2_IM[band] = WFK->cg[kptno][band][pw1][0]*WFK->cg[kptno][band][pw2][1] - WFK->cg[kptno][band][pw1][1]*WFK->cg[kptno][band][pw2][0];
		}
		for(atom_type=0;atom_type<ntypat;atom_type++) {
		  for(l=0;l<l_max_type[atom_type];l++) {
			switch(l) {
			case 0:
			if(i_max_type[atom_type][l]==1) {
			  php = p[pw2][l][0][atom_type-1]*ph[0][0][l][atom_type-1];
			  vg1g2_real_l[atom_type][l] = Plm[pw1][l][0]*php*Plm[pw2][l][0];
			}
			else if(i_max_type[atom_type][l]==2) {
			  php  = p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[1][0][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[0][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[1][1][l][atom_type];
			  vg1g2_real_l[atom_type][l] = Plm[pw1][l][0]*php*Plm[pw2][l][0];
			}
			else if(i_max_type[atom_type][l]==3) {
			  php =  p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[1][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[2][0][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[0][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[1][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[2][1][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[0][2][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[1][2][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[2][2][l][atom_type];
			  vg1g2_real_l[atom_type][l] = Plm[pw1][l][0]*php*Plm[pw2][l][0];
			}
			break;
			case 1:
			if(i_max_type[atom_type][l]==1) {
			  php = p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  vg1g2_real_l[atom_type][l]  = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]);
			}
			else if(i_max_type[atom_type][l]==2) {
			  php =  p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[1][0][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[0][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[1][1][l][atom_type];
			  vg1g2_real_l[atom_type][l] = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]);
			}
			else if(i_max_type[atom_type][l]==3) {
			  php =  p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[1][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[2][0][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[0][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[1][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[2][1][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[0][2][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[1][2][l][atom_type];
			  php += p[pw2][l][2][atom_type]*ph[2][2][l][atom_type];
			  vg1g2_real_l[atom_type][l] = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]);
			}
			break;
			case 2:
			if(i_max_type[atom_type][l]==1) {
			  php = p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  vg1g2_real_l[atom_type][l] = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]+
			  Plm[pw1][l][2]*Plm[pw2][l][2]*cos_mdphi_times_2[2]);
			}
			else if(i_max_type[atom_type][l]==2) {
			  php =  p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  php += p[pw2][l][0][atom_type]*ph[1][0][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[0][1][l][atom_type];
			  php += p[pw2][l][1][atom_type]*ph[1][1][l][atom_type];
			  vg1g2_real_l[atom_type][l] = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]+
			  Plm[pw1][l][2]*Plm[pw2][l][2]*cos_mdphi_times_2[2]);
			}
			break;
			case 3:
			  php = p[pw2][l][0][atom_type]*ph[0][0][l][atom_type];
			  vg1g2_real_l[atom_type][l] = php*(Plm[pw1][l][0]*Plm[pw2][l][0] +
			  Plm[pw1][l][1]*Plm[pw2][l][1]*cos_mdphi_times_2[1]+
			  Plm[pw1][l][2]*Plm[pw2][l][2]*cos_mdphi_times_2[2]+
			  Plm[pw1][l][3]*Plm[pw2][l][3]*cos_mdphi_times_2[3]);
			break;
			}
		  } /* end l loop */
		} /* end atom_type loop*/

		for(atomno=0;atomno<natom;atomno++) {
		  atom_type = ATM->typat[atomno]-1;
		  atom_phase_rad = -(Dga*ATM->xred[0][atomno]+Dgb*ATM->xred[1][atomno]+Dgc*ATM->xred[2][atomno]);
		  table_entry = atom_phase_rad;
		  atom_phase_rad -= table_entry;
		  if (atom_phase_rad<0) atom_phase_rad += 1.000;
		  table_entry = atom_phase_rad/theta_step + 0.5;
		  atomphase_RE = cos_table[table_entry];
		  atomphase_IM = sin_table[table_entry];

		  /*Determine Energy range to store contribution*/
		  for(band=0;band<nband;band++) {
		    bandE = (((WFK->eigen[kptno][band]-fermi)*HATOEV-bandE_min)*Emesh)+0.5;
			dE = floor(bandE);
			if (dE >= nEstep) continue; //skip cont if range too high
			/*Setting Band Occupation and Kpt Weight */
// 		    occ = WFK->occ[kptno][band];
			occ = 1;
 		    wtk = WFK->wtk[kptno];

			/*finding nonlocalE contribution storing by hkl*/
			gkk_real = c1c2_RE[band]*atomphase_RE-c1c2_IM[band]*atomphase_IM;
			c1c2atomphases_real = broad*wtk*occ*gkk_real;
			for (l=0;l<l_max_type[atom_type];l++) {
	          total_nonlocal += 1.0*c1c2atomphases_real*vg1g2_real_l[atom_type][l];
	          ECON->nonlocal[dE] += c1c2atomphases_real*vg1g2_real_l[atom_type][l];
			} //end l->l_max_type loop
		  } //end band->nband loop
		} //end of atom loop

	  } //end pw2 loop
	} //end pw1 loop
	fprintf(flog, "   kpt %d\t nonlocal potential energy = %lf\n", kptno, total_nonlocal);
  } //end of kpt loop
  printf("\tTotal Nonlocal Energy = %lf\n", total_nonlocal);
  fprintf(flog, "Total Nonlocal Energy = %lf\n", total_nonlocal);

  /*free allocated variables*/
  l_max_type = FreeMemory_oneD_int(l_max_type);
  i_max_type = FreeMemory_twoD_int(i_max_type, ntypat);
  ph = FreeMemory_fourD_double(ph, 4, 4, 5);
  c1c2_RE = FreeMemory_oneD_double(c1c2_RE);
  c1c2_IM = FreeMemory_oneD_double(c1c2_IM);
  vg1g2_real_l = FreeMemory_twoD_double(vg1g2_real_l, ntypat);
  gx = FreeMemory_oneD_double(gx);  
  gy = FreeMemory_oneD_double(gy);  
  gz = FreeMemory_oneD_double(gz);  
  ga = FreeMemory_oneD_double(ga);  
  gb = FreeMemory_oneD_double(gb);  
  gc = FreeMemory_oneD_double(gc);  
  g = FreeMemory_oneD_double(g);  
  costheta = FreeMemory_oneD_double(costheta);
  phi = FreeMemory_oneD_double(phi);
  theta = FreeMemory_oneD_double(theta);
  p = FreeMemory_fourD_double(p, PW_max, 4, 4);
  Plm = FreeMemory_threeD_double(Plm, PW_max, 4);

}   //END of Calculate Nonlocal energy contribution function 

