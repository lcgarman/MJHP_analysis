#ifndef psp_functions_H
#define psp_functions_H 

  void xyz2sph(double X, double Y, double Z, double * r, double * theta, double * phi);

  double projp(int l ,int i ,double g1, Pseudopotential* PSP, int typat, double cellvolume);

  void read_PSPdata(char filename[100], Pseudopotential* PSP, AtomicVariables * ATM);

#endif
