#ifndef allocate_memory_H
#define allocate_memory_H  

  void Initialize_NumberGrid(NumberGrid * NG);

  void Initialize_Symmetry(Symmetry * SYM);

  void Initialize_Wavefunction(Wavefunction * WFK);

  void Initialize_BinaryGrid(BinaryGrid * BIN);

  void Initialize_AtomicVariables(AtomicVariables * ATM);

  void Initialize_Pseudopotential(Pseudopotential * PSP);

  void Initialize_VectorIndices(VectorIndices * VECT);

  void Initialize_MottJonesConditions(MottJonesConditions * MJC);

  void Initialize_TwoTheta(TwoTheta * TTH);

  void Initialize_EnergyContribution(EnergyContribution * ECON);
 
  void Initialize_MottJonesPlanewave(MottJonesPlanewave * MJPW);


  int* AllocateMemory_oneD_int(int *array, int dim1);

  int** AllocateMemory_twoD_int(int **array, int dim1, int dim2);

  int*** AllocateMemory_threeD_int(int ***array, int dim1, int dim2, int dim3);

  int**** AllocateMemory_fourD_int(int ****array, int dim1, int dim2, int dim3, int dim4);

  double* AllocateMemory_oneD_double(double *array, int dim1);

  double** AllocateMemory_twoD_double(double **array, int dim1, int dim2);

  double*** AllocateMemory_threeD_double(double ***array, int dim1, int dim2, int dim3);

  double**** AllocateMemory_fourD_double(double ****array, int dim1, int dim2, int dim3, int dim4);

  gsl_complex***** AllocateMemory_fiveD_complex(gsl_complex *****array, int dim1, int dim2, int dim3, int dim4, int dim5);
 
  gsl_complex*** AllocateMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2, int dim3);


  int* FreeMemory_oneD_int(int *array);

  int** FreeMemory_twoD_int(int **array, int dim1);

  int*** FreeMemory_threeD_int(int ***array, int dim1, int dim2);

  int**** FreeMemory_fourD_int(int ****array, int dim1, int dim2, int dim3);

  double* FreeMemory_oneD_double(double *array);

  double** FreeMemory_twoD_double(double **array, int dim1);

  double*** FreeMemory_threeD_double(double ***array, int dim1, int dim2);

  double**** FreeMemory_fourD_double(double ****array, int dim1, int dim2, int dim3);

  gsl_complex*** FreeMemory_threeD_complex(gsl_complex ***array, int dim1, int dim2);

  gsl_complex***** FreeMemory_fiveD_complex(gsl_complex *****array, int dim1, int dim2, int dim3, int dim4);

  void AllocateMemory_Wavefunctions(Wavefunction* WFK); 

  void FreeMemory_Wavefunctions(Wavefunction* WFK); 
  
  void AllocateMemory_PSPvariables(Pseudopotential* PSP, int ntypat);

  void FreeMemory_PSPvariables(Pseudopotential* PSP);

#endif
