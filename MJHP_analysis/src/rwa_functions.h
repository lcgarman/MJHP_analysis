#ifndef rwa_functions_H 
#define rwa_functions_H 

  void read_mjin_header(char filename[100], MottJonesConditions * MJC, FileCabinet* FCAB, FermiSphere* FS);

  void append_mjin(char filename[100], TwoTheta *TTH);

  void print_mjhpHKL_energy(char filename[200], VectorIndices *VECT, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON);

  void read_mjhpHKL_energy(char filename[100], EnergyStep *ESTP, UnitCell *UC, FileCabinet *FCAB, EnergyContribution *ECON);

  void concatinate_hkl_grids(EnergyContribution * ECON, EnergyStep * ESTP); 

  void print_mjhp2theta_energy(char filename[200], TwoTheta *TTH, EnergyStep * ESTP, UnitCell* UC, FileCabinet *FCAB, EnergyContribution * ECON, FermiSphere *FS);

  void print_XSF(char filename[200], UnitCell* UC, NumberGrid* GRD, BinaryGrid* BIN, AtomicVariables* ATM);

  void modify_line(char *line, Modification modifications[], int num_modifications);

  void modify_abinitin(char in_filename[100], char out_filename[100]);

#endif
