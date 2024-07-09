#ifndef hkl_functions_H
#define hkl_functions_H  

  void transform_HKL(MottJonesConditions * MJC); 

  void find_symmetric_hkl(VectorIndices* VECT, Symmetry* SYM, NumberGrid* GRD); 

  void find_MJregion(VectorIndices *VECT, UnitCell *UC); 

  void concatinate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP); 

  void integrate_HKL_potential(EnergyContribution * ECON, EnergyStep * ESTP, AtomicVariables * ATM);

#endif
