#ifndef rwa_functions_H 
#define rwa_functions_H 

  void read_mjin_header(char filename[100], MottJonesConditions * MJC, FileCabinet* FCAB, FermiSphere* FS);

  void print_mjhpHKL_energy(char filename[200], VectorIndices *VECT, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON);

  void print_mjhp2theta_energy(char filename[200], TwoTheta *TTH, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON, FermiSphere *FS);

  void print_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* BIN, AtomicVariables* ATM);

#endif
