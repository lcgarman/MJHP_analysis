#ifndef mjhpHKL_nonlocal_potential_H
#define mjhpHKL_nonlocal_potential_H 

  void mjhpHKL_nonlocal_potential(VectorIndices *VECT, NumberGrid* GRD, EnergyStep* ESTP, Wavefunction* WFK, Pseudopotential* PSP, UnitCell* UC, AtomicVariables* ATM, EnergyContribution *ECON);

  void mjhpREAL_nonlocal_potential(VectorIndices *VECT, NumberGrid* GRD, Wavefunction* WFK, Pseudopotential* PSP, UnitCell* UC, AtomicVariables* ATM, NonlocalEnergyReal* NLER);
  
#endif
