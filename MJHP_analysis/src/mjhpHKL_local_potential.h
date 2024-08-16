#ifndef mjhpHKL_local_potential_H 
#define mjhpHKL_local_potential_H 

  void mjhpHKL_local_potential(VectorIndices *VECT, NumberGrid *GRD, EnergyStep *ESTP, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, EnergyContribution *ECON);

  void mjhpREAL_local_potential(VectorIndices *VECT, NumberGrid *GRD, Wavefunction *WFK, UnitCell *UC, BinaryGrid *BIN, BinaryGrid *LOCAL);

#endif
