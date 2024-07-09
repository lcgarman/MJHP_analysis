#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "globals.h"
#include "structures.h"
#include "allocate_memory.h"
#include "psp_functions.h"

void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi)
{
  /* Z = r cos(theta), X = r sin(theta) cos(phi) Y = r sin(theta)sin(phi) */
  *r = pow((pow(X,2.0)+pow(Y,2.0)+pow(Z,2.0)),0.5);
  if(*r>0.0) {
	*theta=acos(Z/ *r);
	if(fabs(sin(*theta))>0.0) {
	  *phi=acos(X/(*r * sin(*theta)));
	  if(X/(*r * sin(*theta))>1.0) {
		*phi=acos(1.0);
	  }
	  if(X/(*r * sin(*theta))<-1.0) {
		*phi=acos(-1.0);
	  }
	  if(Y<0.0) *phi=-*phi;
	}
	else *phi=0;
  }
  else {
	*phi=0.0;
	*theta=0.0;
  }
}//END of xyz2sph function
        
double projp(int l ,int i ,double g1, Pseudopotential* PSP, int typat, double cellvolume)
{
  double p=0.0;
  typat--;
  if((l==0)&&(i==0)) {
	p=4.0*pow(2.0*pow(PSP->rrs[typat],3.0),0.5)*PI_5_4/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrs[typat],2.0)));
  }
  if((l==0)&&(i==1)) {
	p=8.0*pow(2.0*pow(PSP->rrs[typat],3.0)/15.0,0.5)*PI_5_4*(3.0-pow(g1*PSP->rrs[typat],2.0))/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrs[typat],2.0)));	
  }
  if((l==0)&&(i==2)) {
	p=16.0*pow(2.0*pow(PSP->rrs[typat],3.0)/105.0,0.5)*PI_5_4*(15.0-10.0*pow(g1*PSP->rrs[typat],2.0)+pow(g1*PSP->rrs[typat],4.0))/(3.0*pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrs[typat],2.0)));	
  }
  if((l==1)&&(i==0)) {
	p=8.0*pow(pow(PSP->rrp[typat],5.0)/3.0,0.5)*PI_5_4*(g1)/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrp[typat],2.0)));	
  }
  if((l==1)&&(i==1)) {
	p=16.0*pow(pow(PSP->rrp[typat],5.0)/105.0,0.5)*PI_5_4*(g1)*(5.0-pow(g1*PSP->rrp[typat],2.0))/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrp[typat],2.0)));	
  }
  if((l==1)&&(i==2)) {
	p=32.0*pow(pow(PSP->rrp[typat],5.0)/1155.0,0.5)*PI_5_4*(g1)*(35.0-14.0*pow(g1*PSP->rrp[typat],2.0)+pow(g1*PSP->rrp[typat],4.0))/(3.0*pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrp[typat],2.0)));	
  }
  if((l==2)&&(i==0)) {
	p=8.0*pow(2*pow(PSP->rrd[typat],7.0)/15.0,0.5)*PI_5_4*(g1*g1)/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrd[typat],2.0)));	
  }
  if((l==2)&&(i==1)) {
	p=16.0*pow(2*pow(PSP->rrd[typat],7.0)/105.0,0.5)*PI_5_4*(g1*g1)*(7.0-pow(g1*PSP->rrd[typat],2.0))/(3.0*pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrd[typat],2.0)));	
  }	
  if((l==3)&&(i==0)) {
	p=16.0*pow(pow(PSP->rrf[typat],9.0)/105.0,0.5)*PI_5_4*(g1*g1*g1)/(pow(cellvolume,0.5)*exp(0.5*pow(g1*PSP->rrf[typat],2.0)));
  }
  // if(p<0.0) printf("YIKES!");
  return(p);
}

void read_PSPdata(char filename[100], Pseudopotential* PSP, AtomicVariables * ATM)
{
  FILE* faout;
  int stop;
  int stop_ntypat;
  int check; 
  int ntypat;
  int typat;
  char str [100];
  double rloc;
  double cc1, cc2, cc3, cc4;
  double rrs, rrp, rrd, rrf;
  double k11p, k22p, k33p;
  double k11d, k22d, k33d;
  double k11f, k22f, k33f;
  double h11s, h22s, h33s;
  double h11p, h22p, h33p;
  double h11d, h22d, h33d;
  double h11f, h22f, h33f;
  
  ntypat = ATM->ntypat;

  printf( "\nReading %s for PSP data\n", filename);
  /*open abinit *.out file in read mode*/
  faout = fopen(filename, "r");
  if(faout==NULL) {
    printf("%s not found. \n", filename);
    exit(0);
  }  

  /*Allocate Memory for psp variables*/
  AllocateMemory_PSPvariables(PSP, ntypat);
  /*end of allocation*/

  stop=0;
  typat = 0;
  /*Begin reading in from the out file*/
  /*Store in h[#][#][l][atomtype]*/
  while(stop==0) {
	check=fscanf(faout,"%s",str);
	if(check==EOF) {
	  stop=1;
	}
	if(strcmp(str,"pspini:")==0) {
	  stop_ntypat=0;
	  while(stop_ntypat==0) {
		check=fscanf(faout,"%s",str);
		if(strcmp(str,"rloc=")==0) {
		  fscanf(faout, "%lf", &rloc);
		  PSP->rloc[typat] = rloc;
		}
		if(strcmp(str,"cc1")==0) {
		  fscanf(faout, " =  %lf", &cc1);
		  PSP->cc1[typat] = cc1;
		}
		if(strcmp(str,"cc2")==0) {
		  fscanf(faout, " =  %lf", &cc2);
		  PSP->cc2[typat] = cc2;
		}
		if(strcmp(str,"cc3")==0) {
		  fscanf(faout, " =  %lf", &cc3);
		  PSP->cc3[typat] = cc3;
		}
		if(strcmp(str,"cc4")==0) {
		  fscanf(faout, " =  %lf", &cc4);
		  PSP->cc4[typat] = cc4;
		}
		if(strcmp(str,"rrs")==0) {
		  fscanf(faout, " =  %lf", &rrs);
		  PSP->rrs[typat] = rrs;
		}
		if(strcmp(str,"h11s=")==0) {
		  fscanf(faout, "%lf", &h11s);
		  PSP->h[0][0][0][typat] = h11s;
		} 
        if(strcmp(str, "h22s=")==0) {
          fscanf(faout, "%lf", &h22s);
		  PSP->h[1][1][0][typat] = h22s;
        } 
        if(strcmp(str, "h33s=")==0) {
          fscanf(faout, "%lf", &h33s);
		  PSP->h[2][2][0][typat] = h33s;
        } 
		if(strcmp(str,"rrp")==0) {
		  fscanf(faout, " =  %lf", &rrp);
          PSP->rrp[typat] = rrp;
        } 
        if(strcmp(str, "h11p=")==0) {
          fscanf(faout, "%lf", &h11p);
		  PSP->h[0][0][1][typat] = h11p;
        } 
        if(strcmp(str, "h22p=")==0) {
          fscanf(faout, "%lf", &h22p);
		  PSP->h[1][1][1][typat] = h22p;
        } 
        if(strcmp(str, "h33p=")==0) {
          fscanf(faout, "%lf", &h33p);
		  PSP->h[2][2][1][typat] = h33p;
        } 
        if(strcmp(str, "k11p=")==0) {
          fscanf(faout, "%lf", &k11p);
          PSP->k11p[typat] = k11p;
        } 
        if(strcmp(str, "k22p=")==0) {
          fscanf(faout, "%lf", &k22p);
          PSP->k22p[typat] = k22p;
        } 
        if(strcmp(str, "k33p=")==0) {
          fscanf(faout, "%lf", &k33p);
          PSP->k33p[typat] = k33p;
        } 
		if(strcmp(str,"rrd")==0) {
		  fscanf(faout, " =  %lf", &rrd);
          PSP->rrd[typat] = rrd;
        } 
        if(strcmp(str, "h11d=")==0) {
          fscanf(faout, "%lf", &h11d);
		  PSP->h[0][0][2][typat] = h11d;
        } 
        if(strcmp(str, "h22d=")==0) {
          fscanf(faout, "%lf", &h22d);
		  PSP->h[1][1][2][typat] = h22d;
        } 
        if(strcmp(str, "h33d=")==0) {
          fscanf(faout, "%lf", &h33d);
		  PSP->h[2][2][2][typat] = h33d;
        } 
        if(strcmp(str, "k11d=")==0) {
          fscanf(faout, "%lf", &k11d);
          PSP->k11d[typat] = k11d;
        } 
        if(strcmp(str, "k22d=")==0) {
          fscanf(faout, "%lf", &k22d);
          PSP->k22d[typat] = k22d;
        } 
        if(strcmp(str, "k33d=")==0) {
          fscanf(faout, "%lf", &k33d);
          PSP->k33d[typat] = k33d;
        } 
		if(strcmp(str,"rrf")==0) {
		  fscanf(faout, " =  %lf", &rrf);
          PSP->rrf[typat] = rrf;
        } 
        if(strcmp(str, "h11f=")==0) {
          fscanf(faout, "%lf", &h11f);
		  PSP->h[0][0][3][typat] = h11f;
        } 
        if(strcmp(str, "h22f=")==0) {
          fscanf(faout, "%lf", &h22f);
		  PSP->h[1][1][3][typat] = h22f;
        } 
        if(strcmp(str, "h33f=")==0) {
          fscanf(faout, "%lf", &h33f);
		  PSP->h[2][2][3][typat] = h33f;
        } 
        if(strcmp(str, "k11f=")==0) {
          fscanf(faout, "%lf", &k11f);
          PSP->k11f[typat] = k11f;
        } 
        if(strcmp(str, "k22f=")==0) {
          fscanf(faout, "%lf", &k22f);
          PSP->k22f[typat] = k22f;
        } 
        if(strcmp(str, "k33f=")==0) {
          fscanf(faout, "%lf", &k33f);
          PSP->k33f[typat] = k33f;
        } 
        if(strcmp(str,"COMMENT")==0) {
		  typat++;
		  if (typat==ntypat) stop=1;
		  stop_ntypat = 1;
        }

	  } //END while (stop_ntypat==0) loop
	} //END search for pspini
  } //END while (stop==0) loop

  printf( "Number of Atom Types = %d\n", ntypat);
  ATM->ntypat = ntypat;
     /*Now define symmetric off diagnol elements in h*/
  for (typat=0;typat<ntypat;typat++) {
	PSP->h[0][1][0][typat]=-0.5*pow(3.0/5.0,0.5)*PSP->h[1][1][0][typat];
	PSP->h[1][0][0][typat]=PSP->h[0][1][0][typat];
	PSP->h[0][2][0][typat] = (0.5)*pow(5.0/21.0,0.5)*PSP->h[2][2][0][typat];
	PSP->h[2][0][0][typat] = PSP->h[0][2][0][typat];
	PSP->h[1][2][0][typat] = (-0.5)*pow(100.0/63.0,0.5)*PSP->h[2][2][0][typat];
	PSP->h[2][1][0][typat] = PSP->h[1][2][0][typat];
	PSP->h[0][1][1][typat] = (-0.5)*pow(5.0/7.0,0.5)*PSP->h[1][1][1][typat];
	PSP->h[1][0][1][typat] = PSP->h[0][1][1][typat];
	PSP->h[0][2][1][typat] = (1.0/6.0)*pow(35.0/11.0,0.5)*PSP->h[2][2][1][typat];
	PSP->h[2][0][1][typat] = PSP->h[0][2][1][typat];
	PSP->h[1][2][1][typat] = (-1.0/6.0)*(14.0/pow(11.0, 0.5))*PSP->h[2][2][1][typat];
	PSP->h[2][1][1][typat] = PSP->h[1][2][1][typat];
	PSP->h[0][1][2][typat] = (-0.5)*pow(7.0/9.0, 0.5)*PSP->h[1][1][2][typat];
	PSP->h[1][0][2][typat] = PSP->h[0][1][2][typat];
	PSP->h[0][2][2][typat] = (0.5)*pow(63.0/143.0, 0.5)*PSP->h[2][2][2][typat];
	PSP->h[2][0][2][typat] = PSP->h[0][2][2][typat];
	PSP->h[1][2][2][typat] = (-0.5)*(18.0/pow(143.0,0.5))*PSP->h[2][2][2][typat]; 
	PSP->h[2][1][2][typat] = PSP->h[1][2][2][typat];
	
	printf( "Pseudopotential for atomtype # %d\n", typat);
	printf( "   rloc= %lf \n",PSP->rloc[typat]);
	printf( "   cc1 = %lf; cc2 = %lf; cc3 = %lf; cc4 = %lf \n",PSP->cc1[typat],PSP->cc2[typat],PSP->cc3[typat],PSP->cc4[typat]);
	printf( "   rrs = %lf; h11s= %lf; h22s= %lf; h33s= %lf \n",PSP->rrs[typat],PSP->h[0][0][0][typat],PSP->h[1][1][0][typat],PSP->h[2][2][0][typat]);
	printf( "   rrp = %lf; h11p= %lf; h22p= %lf; h33p= %lf \n",PSP->rrp[typat],PSP->h[0][0][1][typat],PSP->h[1][1][1][typat],PSP->h[2][2][1][typat]);
	printf( "   rrp = %lf; k11p= %lf; k22p= %lf; k33p= %lf \n",PSP->rrp[typat],PSP->k11p[typat],PSP->k22p[typat],PSP->k33p[typat]);
	printf( "   rrd = %lf; h11d= %lf; h22d= %lf; h33d= %lf \n",PSP->rrd[typat],PSP->h[0][0][2][typat],PSP->h[1][1][2][typat],PSP->h[2][2][2][typat]);
	printf( "   rrd = %lf; k11d= %lf; k22d= %lf; k33d= %lf \n",PSP->rrd[typat],PSP->k11d[typat],PSP->k22d[typat],PSP->k33d[typat]);
	printf( "   rrf = %lf; h11f= %lf; h22f= %lf; h33f= %lf \n",PSP->rrf[typat],PSP->h[0][0][3][typat],PSP->h[1][1][3][typat],PSP->h[2][2][3][typat]);
	printf( "   rrf = %lf; k11f= %lf; k22f= %lf; k33f= %lf \n",PSP->rrf[typat],PSP->k11f[typat],PSP->k22f[typat],PSP->k33f[typat]);
	printf( "\n");
	/*	printf("h000 = %lf\n", PSP->h[0][0][0][typat]);
	printf("h100 = %lf\n", PSP->h[1][0][0][typat]);
	printf("h010 = %lf\n", PSP->h[0][1][0][typat]);
	printf("h110 = %lf\n", PSP->h[1][1][0][typat]);
	printf("h200 = %lf\n", PSP->h[2][0][0][typat]);
	printf("h020 = %lf\n", PSP->h[0][2][0][typat]);
	printf("h220 = %lf\n", PSP->h[2][2][0][typat]);
	printf("h120 = %lf\n", PSP->h[1][2][0][typat]);
	printf("h210 = %lf\n", PSP->h[2][1][0][typat]);
	printf("h001 = %lf\n", PSP->h[0][0][1][typat]);
	printf("h101 = %lf\n", PSP->h[1][0][1][typat]);
	printf("h011 = %lf\n", PSP->h[0][1][1][typat]);
	printf("h111 = %lf\n", PSP->h[1][1][1][typat]);
	printf("h201 = %lf\n", PSP->h[2][0][1][typat]);
	printf("h021 = %lf\n", PSP->h[0][2][1][typat]);
	printf("h221 = %lf\n", PSP->h[2][2][1][typat]);
	printf("h121 = %lf\n", PSP->h[1][2][1][typat]);
	printf("h211 = %lf\n", PSP->h[2][1][1][typat]);
	printf("h002 = %lf\n", PSP->h[0][0][2][typat]);
	printf("h102 = %lf\n", PSP->h[1][0][2][typat]);
	printf("h012 = %lf\n", PSP->h[0][1][2][typat]);
	printf("h112 = %lf\n", PSP->h[1][1][2][typat]);
	printf("h202 = %lf\n", PSP->h[2][0][2][typat]);
	printf("h022 = %lf\n", PSP->h[0][2][2][typat]);
	printf("h222 = %lf\n", PSP->h[2][2][2][typat]);
	printf("h122 = %lf\n", PSP->h[1][2][2][typat]);
	printf("h212 = %lf\n", PSP->h[2][1][2][typat]);
*/  }
  fclose(faout);

} //END of Read_PSPinfo function
  
