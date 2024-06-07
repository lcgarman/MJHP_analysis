#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "structures.h"
#include "allocate_memory.h"

void Initialize_NumberGrid(NumberGrid * NG) 
{
  NG->h_grid = NULL;
  NG->k_grid = NULL;
  NG->l_grid = NULL;
}

void Initialize_Symmetry(Symmetry * SYM) 
{
  SYM->symrel = NULL;
  SYM->mult = NULL;
  SYM->tnons = NULL;
  SYM->symor = NULL;
}
 
void Initialize_Wavefunction(Wavefunction * WFK)
{
  WFK->npw = NULL;
  WFK->kpt = NULL;
  WFK->wtk = NULL;
  WFK->kg = NULL;
  WFK->eigen = NULL;
  WFK->occ = NULL;
  WFK->cg = NULL;
}

void Initialize_BinaryGrid(BinaryGrid * BIN)
{
  BIN->real_grid = NULL;
  BIN->rec_grid = NULL;
  BIN->cc_rec_grid = NULL;
}

void Initialize_AtomicVariables(AtomicVariables * ATM) 
{
  ATM->xred = NULL;
  ATM->typat = NULL;
  ATM->atomicno = NULL;
  ATM->Xcart = NULL;
  ATM->Ycart = NULL;
  ATM->Zcart = NULL;
}

void Initialize_Pseudopotential(Pseudopotential * PSP)
{
  PSP->rloc = NULL;
  PSP->cc1 = NULL;
  PSP->cc2 = NULL;
  PSP->cc3 = NULL;
  PSP->cc4 = NULL;
  PSP->rrs = NULL;
  PSP->rrp = NULL;
  PSP->rrd = NULL;
  PSP->rrf = NULL;
  PSP->k11p = NULL;
  PSP->k22p = NULL;
  PSP->k33p = NULL;
  PSP->k11d = NULL;
  PSP->k22d = NULL;
  PSP->k33d = NULL;
  PSP->k11f = NULL;
  PSP->k22f = NULL;
  PSP->k33f = NULL;
  PSP->h = NULL;
}

void Initialize_VectorIndices(VectorIndices * VECT) 
{
  VECT->H_arr = NULL;
  VECT->K_arr = NULL;
  VECT->L_arr = NULL;
  VECT->H_posarr = NULL;
  VECT->K_posarr = NULL;
  VECT->L_posarr = NULL;
}

void Initialize_MottJonesConditions(MottJonesConditions * MJC)
{
  MJC->jzH = NULL;
  MJC->jzK = NULL;
  MJC->jzL = NULL;
  MJC->JZkpt = NULL;
  MJC->JZk_sym = NULL;
  MJC->JZmult = NULL;
  MJC->JZwtk = NULL;
  MJC->BZk = NULL;
  MJC->Hpw = NULL;
  MJC->Kpw = NULL;
  MJC->Lpw = NULL;
}

void Initialize_TwoTheta(TwoTheta * TTH)
{
  TTH->two_theta = NULL;
  TTH->d_hkl = NULL;
  TTH->F_hkl = NULL;
  TTH->rflc_H = NULL;
  TTH->rflc_K = NULL;
  TTH->rflc_L = NULL;
  TTH->hpw = NULL;
  TTH->kpw = NULL;
  TTH->lpw = NULL;
  TTH->rflc_mult = NULL;
  TTH->kx = NULL;
  TTH->ky = NULL;
  TTH->kz = NULL;
  TTH->BZkpt = NULL;
  TTH->BZkpt_sym = NULL;
  TTH->hpw_sym = NULL;
  TTH->kpw_sym = NULL;
  TTH->lpw_sym = NULL;
  TTH->rflc_H_sym = NULL;
  TTH->rflc_K_sym = NULL;
  TTH->rflc_L_sym = NULL;
}

void Initialize_EnergyContribution(EnergyContribution * ECON)
{
  ECON->local = NULL;
  ECON->nonlocal = NULL;
  ECON->total_potential = NULL;
  ECON->rflc_local = NULL;
  ECON->rflc_nonlocal = NULL;
  ECON->rflc_total = NULL;
}


int* AllocateMemory_oneD_int(int *array, int dim1)
{
  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }
    
  /*allocate 1D array to pointer*/
  array = malloc(dim1 * sizeof(int));
  /*check allocation was successful*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*return the array at end of function*/
  return array;
} 

int** AllocateMemory_twoD_int(int **array, int dim1, int dim2)
{
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int*));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed (1D/2)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed(2D/2)\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  /*return the array at end of function*/
  return array;
} 

int*** AllocateMemory_threeD_int(int ***array, int dim1, int dim2, int dim3)
{
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(int));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

int**** AllocateMemory_fourD_int(int ****array, int dim1, int dim2, int dim3, int dim4)
{
  int i, j, k, l;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(int***));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(int**));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(int*));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }

      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(int));
        /*check allocation*/
		if (array[i][j][k] == NULL) {
		  printf("ERROR: Memory Allocation Failed\n");
		  for (l=0;l<k;l++) {
			free(array[i][j][l]);
		  }
		  free(array[i][j]);
		  free(array[i]);
		  free(array);
		  exit(0);
		}
	  }
    }
  }

  /*return the array at end of function*/
  return array;
} 

double* AllocateMemory_oneD_double(double *array, int dim1)
{
  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
printf("NOT NULL\n");
  }

  /*allocate 1D array to pointer*/
  array = malloc(dim1 * sizeof(double));
  /*check allocation was successful*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*return the array at end of function*/
  return array;
} 

double** AllocateMemory_twoD_double(double **array, int dim1, int dim2)
{
  int i, j;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of array*/
  array = malloc(dim1 * sizeof(double*));
  /*check allocation for every dimension*/
  if (array==NULL) {
    printf("ERROR: Memory Allocation Failed (1D)\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }
  }

  return array;
}

double*** AllocateMemory_threeD_double(double ***array, int dim1, int dim2, int dim3)
{
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 

double**** AllocateMemory_fourD_double(double ****array, int dim1, int dim2, int dim3, int dim4)
{
  int i, j, k, l;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(double***));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(double**));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(double*));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
      
      /*allocate fourth dimension*/
      for (k=0;k<dim3;k++) {
        array[i][j][k] = malloc(dim4 * sizeof(double));
        /*check allocation*/
        if (array[i][j][k] == NULL) {
          printf("ERROR: Memory Allocation Failed\n");
          for (l=0;l<k;l++) {
            free(array[i][j][l]);
          }
          free(array[i][j]);
          free(array[i]);
          free(array);
          exit(0);
        }
	  }
    }
  }

  /*return the array at end of function*/
  return array;
} 

gsl_complex*** AllocateMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2, int dim3)
{
  int i, j, k;

  /*if array is already allocated leave this function*/
  if (array != NULL) {
    return array;
  }

  /*allocate first dimension of 3D array to pointer*/
  array = malloc(dim1 * sizeof(gsl_complex**));
  /*check if allocation was successful*/
  if (array == NULL) {
    printf("ERROR: Memory Allocation Failed\n");
    exit(0);
  }

  /*allocate second dimension*/ 
  for (i=0;i<dim1;i++) {
    array[i] = malloc(dim2 * sizeof(gsl_complex*));
    /*check allocation*/
    if (array[i] == NULL) {
      printf("ERROR: Memory Allocation Failed\n");
      /*free previously allocated memory*/
      for (j=0;j<i;j++) {
        free(array[j]);
      }
      free(array);
      exit(0);
    }

    /*allocate third dimension*/
    for (j=0;j<dim2;j++) {
      array[i][j] = malloc(dim3 * sizeof(gsl_complex));
      /*check allocation*/
      if (array[i][j] == NULL) {
        printf("ERROR: Memory Allocation Failed\n");
        for (k=0;k<j;k++) {
          free(array[i][k]);
        }
        free(array[i]);
        free(array);
		exit(0);
      }
    }
  }

  /*return the array at end of function*/
  return array;
} 



int* FreeMemory_oneD_int(int *array)
{
  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  free(array);
  array = NULL;
  return array;
}

int** FreeMemory_twoD_int(int **array, int dim1)
{
  int i;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

int*** FreeMemory_threeD_int(int ***array, int dim1, int dim2)
{
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

int**** FreeMemory_fourD_int(int ****array, int dim1, int dim2, int dim3)
{
  int i, j, k;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

double* FreeMemory_oneD_double(double *array)
{
  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  free(array);
  array = NULL;
  return array;
}

double** FreeMemory_twoD_double(double **array, int dim1)
{
  int i;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

double*** FreeMemory_threeD_double(double ***array, int dim1, int dim2)
{
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;

  return array;
}

double**** FreeMemory_fourD_double(double ****array, int dim1, int dim2, int dim3)
{
  int i, j, k;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }

  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      for (k=0;k<dim3;k++) {
        free(array[i][j][k]);
      }
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}

gsl_complex*** FreeMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2)
{
  int i, j;

  /*if array is NULL leave this function*/
  if (array == NULL) {
    return array;
  }
  
  for (i=0;i<dim1;i++) {
    for (j=0;j<dim2;j++) {
      free(array[i][j]);
    }
    free(array[i]);
  }
  free(array);
  array = NULL;
  return array;
}


void AllocateMemory_Wavefunctions(Wavefunction* WFK) 
{
  int nkpt;
  int nband;
  int npw;
  int i, j, k;

  nkpt = WFK->nkpt;
  nband = WFK->nband;

  if (WFK->cg != NULL) {
    printf("ERROR: wavefunction coefficients not empty before allocation\n");
  }
  if (WFK->kg != NULL) {
    printf("ERROR: wavefunction coefficients not empty before allocation\n");
  }

  /*Allocating Memory for Wavefunction Coefficients cg*/
  WFK->cg = malloc(nkpt * sizeof(double***));
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    WFK->cg[i] = malloc(nband * sizeof(double**));
    for (j=0;j<nband;j++) {
      WFK->cg[i][j] = malloc(npw * sizeof(double*));
	  for (k=0;k<npw;k++) {
        WFK->cg[i][j][k] = malloc(2 * sizeof(double));
	  }
	}
  }
  /*Checking if Memory was allocated successfully*/
  if (WFK->cg == NULL) {
	printf("ERROR: Memory Allocation Failed for Wavefunction Coefficients (cg) \n");
    exit(0);
  } 

  /*Allocating Memory for Reduced plane wave Coordinates (kg)*/
  WFK->kg = malloc(nkpt * sizeof(int**));
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    WFK->kg[i] = malloc(npw * sizeof(int*));
    for (k=0;k<npw;k++) {
      WFK->kg[i][k] = malloc(3 * sizeof(int));
	}
  }
  if (WFK->kg == NULL) {
	printf("\nERROR: Memory Allocation Failed for PW Coordinates (kg) \n");
    exit(0);
  }
  
}
void FreeMemory_Wavefunctions(Wavefunction* WFK) 
{
  int nkpt;
  int nband;
  int npw;
  int i, j, k;

  nkpt = WFK->nkpt;
  nband = WFK->nband;

  if ((WFK->kg == NULL)) {
    return;
  }
  if ((WFK->cg == NULL)) {
    return;
  }

  /*free planewave reduced coordinates*/
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (k=0;k<npw;k++) {
      free(WFK->kg[i][k]);
    }
    free(WFK->kg[i]);
  }
  free(WFK->kg);
  WFK->kg = NULL;

  /*free wavefunction coefficients*/
  for (i=0;i<nkpt;i++) {
    npw = WFK->npw[i];
    for (j=0;j<nband;j++) {
      for (k=0;k<npw;k++) {
        free(WFK->cg[i][j][k]);
      }
      free(WFK->cg[i][j]);
    }
    free(WFK->cg[i]);
  }
  free(WFK->cg);
  WFK->cg = NULL;

}

void AllocateMemory_PSPvariables(Pseudopotential* PSP, int ntypat)
{
  /*Initialize psp variables*/
  PSP->rloc = NULL; 
  PSP->cc1 = NULL; 
  PSP->cc2 = NULL; 
  PSP->cc3 = NULL; 
  PSP->cc4 = NULL; 
  PSP->rrs = NULL; 
  PSP->rrp = NULL; 
  PSP->rrd = NULL; 
  PSP->rrf = NULL; 
  PSP->k11d = NULL; 
  PSP->k22d = NULL; 
  PSP->k33d = NULL; 
  PSP->k11f = NULL; 
  PSP->k22f = NULL; 
  PSP->k33f = NULL; 
  PSP->h = NULL; 

  /*allocate Memory*/
  PSP->rloc = AllocateMemory_oneD_double(PSP->rloc, ntypat); 
  PSP->cc1 = AllocateMemory_oneD_double(PSP->cc1, ntypat); 
  PSP->cc2 = AllocateMemory_oneD_double(PSP->cc2, ntypat); 
  PSP->cc3 = AllocateMemory_oneD_double(PSP->cc3, ntypat); 
  PSP->cc4 = AllocateMemory_oneD_double(PSP->cc4, ntypat); 
  PSP->rrs = AllocateMemory_oneD_double(PSP->rrs, ntypat); 
  PSP->rrp = AllocateMemory_oneD_double(PSP->rrp, ntypat); 
  PSP->rrd = AllocateMemory_oneD_double(PSP->rrd, ntypat); 
  PSP->rrf = AllocateMemory_oneD_double(PSP->rrf, ntypat); 
  PSP->k11p = AllocateMemory_oneD_double(PSP->k11p, ntypat); 
  PSP->k22p = AllocateMemory_oneD_double(PSP->k22p, ntypat); 
  PSP->k33p = AllocateMemory_oneD_double(PSP->k33p, ntypat); 
  PSP->k11d = AllocateMemory_oneD_double(PSP->k11d, ntypat); 
  PSP->k22d = AllocateMemory_oneD_double(PSP->k22d, ntypat); 
  PSP->k33d = AllocateMemory_oneD_double(PSP->k33d, ntypat); 
  PSP->k11f = AllocateMemory_oneD_double(PSP->k11f, ntypat); 
  PSP->k22f = AllocateMemory_oneD_double(PSP->k22f, ntypat); 
  PSP->k33f = AllocateMemory_oneD_double(PSP->k33f, ntypat); 
  PSP->h = AllocateMemory_fourD_double(PSP->h, 4, 4, 5, ntypat);
} 
void FreeMemory_PSPvariables(Pseudopotential* PSP)
{
  PSP->rloc = FreeMemory_oneD_double(PSP->rloc); 
  PSP->cc1 = FreeMemory_oneD_double(PSP->cc1); 
  PSP->cc2 = FreeMemory_oneD_double(PSP->cc2); 
  PSP->cc3 = FreeMemory_oneD_double(PSP->cc3); 
  PSP->cc4 = FreeMemory_oneD_double(PSP->cc4); 
  PSP->rrs = FreeMemory_oneD_double(PSP->rrs); 
  PSP->rrp = FreeMemory_oneD_double(PSP->rrp); 
  PSP->rrd = FreeMemory_oneD_double(PSP->rrd); 
  PSP->rrf = FreeMemory_oneD_double(PSP->rrf); 
  PSP->k11p = FreeMemory_oneD_double(PSP->k11p); 
  PSP->k22p = FreeMemory_oneD_double(PSP->k22p); 
  PSP->k33p = FreeMemory_oneD_double(PSP->k33p); 
  PSP->k11d = FreeMemory_oneD_double(PSP->k11d); 
  PSP->k22d = FreeMemory_oneD_double(PSP->k22d); 
  PSP->k33d = FreeMemory_oneD_double(PSP->k33d); 
  PSP->k11f = FreeMemory_oneD_double(PSP->k11f); 
  PSP->k22f = FreeMemory_oneD_double(PSP->k22f); 
  PSP->k33f = FreeMemory_oneD_double(PSP->k33f); 
  PSP->h = FreeMemory_fourD_double(PSP->h, 4, 4, 5);
}
 
